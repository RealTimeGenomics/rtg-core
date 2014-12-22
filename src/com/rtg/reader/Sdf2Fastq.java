/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.reader;


import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.launcher.CommonFlags.STDIO_NAME;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.Locale;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.mode.DnaUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;

/**
 * This class takes format directory and converts to FASTQ format
 */
public final class Sdf2Fastq extends AbstractCli {

  static final String MODULE_NAME = "sdf2fastq";

  static final String INPUT = "input";
  static final String OUTPUT_FILE = "output";
  static final String LINE_FLAG = "line-length";
  static final String RENAME = "Xrename";
  static final String DEFAULT_QUALITY = "default-quality";

  static final String START_SEQUENCE = "start-id";
  static final String END_SEQUENCE = "end-id";

  static final String FASTQ_EXT_LGE = ".fastq";
  static final String FASTQ_EXT_SRT = ".fq";

  private static String sScore = null;

  private static final Validator VALIDATOR = new Validator() {
    @Override
    public boolean isValid(final CFlags flags) {
      if ((Integer) flags.getValue(LINE_FLAG) < 0) {
        Diagnostic.error(ErrorType.EXPECTED_NONNEGATIVE, LINE_FLAG);
        return false;
      }
      final File input = (File) flags.getValue(INPUT);
      if (!CommonFlags.validateSDF(input)) {
        return false;
      }
      if (flags.isSet(DEFAULT_QUALITY)) {
        final int qual = (Integer) flags.getValue(DEFAULT_QUALITY);
        if (qual < 0) {
          Diagnostic.error(ErrorType.INVALID_MIN_INTEGER_FLAG_VALUE, DEFAULT_QUALITY, Integer.toString(qual), "0");
          return false;
        }
        if (qual > 63) {
          Diagnostic.error(ErrorType.INVALID_MAX_INTEGER_FLAG_VALUE, DEFAULT_QUALITY, Integer.toString(qual), "63");
          return false;
        }
      }
      return CommonFlags.validateStartEnd(flags, START_SEQUENCE, END_SEQUENCE);
    }
  };

  static void warnIfDeletionFails(final File sequenceOutputFile) {
    if (sequenceOutputFile != null && sequenceOutputFile.isFile()) {
      if (sequenceOutputFile.exists() && !sequenceOutputFile.delete()) {
        Diagnostic.warning(WarningType.INFO_WARNING, "Could not delete file \"" + sequenceOutputFile.getAbsolutePath() + "\"");
      }
    }
  }

  /**
   * @return current name of the module
   */
  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  /**
   * Convert an SDF to FASTQ
   *
   * @param out where the FASTQ file goes
   * @param err where error messages go
   * @throws IOException should an IO Error occur
   * @return 0 for success, 1 for failure
   */
  @Override
  protected int mainExec(final OutputStream out, final PrintStream err) throws IOException {
    final File preReadDir = (File) mFlags.getValue(INPUT);
    String output = mFlags.getValue(OUTPUT_FILE).toString();
    String ext = FASTQ_EXT_LGE;
    if (FileUtils.isGzipFilename(output)) {
      output = output.substring(0, output.length() - FileUtils.GZ_SUFFIX.length());
    }
    if (output.toLowerCase(Locale.getDefault()).endsWith(FASTQ_EXT_LGE)) {
      ext = output.substring(output.length() - FASTQ_EXT_LGE.length());
      output = output.substring(0, output.length() - FASTQ_EXT_LGE.length());
    } else if (output.toLowerCase(Locale.getDefault()).endsWith(FASTQ_EXT_SRT)) {
      ext = output.substring(output.length() - FASTQ_EXT_SRT.length());
      output = output.substring(0, output.length() - FASTQ_EXT_SRT.length());
    }
    final File outputFile = new File(output);
    final boolean isPaired = ReaderUtils.isPairedEndDirectory(preReadDir);
    try {
      if (isPaired) {
        if (CommonFlags.isStdio(output)) {
          err.println("Sending paired-end data to stdout is not supported.");
          return 1;
        }
        final LongRange calculatedRegion = doPrereadToFastq(ReaderUtils.getLeftEnd(preReadDir), outputFile.getAbsolutePath() + "_1" + ext, out, null);
        doPrereadToFastq(ReaderUtils.getRightEnd(preReadDir), outputFile.getAbsolutePath() + "_2" + ext, out, calculatedRegion);
      } else {
        doPrereadToFastq(preReadDir, CommonFlags.isStdio(output) ? STDIO_NAME : (outputFile.getAbsolutePath() + ext), out, null);
      }
    } catch (InvalidParamsException e) {
      e.printErrorNoLog();
      return 1;
    }
    return 0;
  }

  //Calculated region for sloppy end
  private LongRange doPrereadToFastq(File preReadDir, String output, OutputStream out, LongRange calculatedRegion) throws IOException, InvalidParamsException {
    final long startId = mFlags.isSet(START_SEQUENCE) ? (Long) mFlags.getValue(START_SEQUENCE) : LongRange.MISSING;
    final long endId = calculatedRegion != null ? calculatedRegion.getEnd() : (mFlags.isSet(END_SEQUENCE) ? (Long) mFlags.getValue(END_SEQUENCE) : LongRange.MISSING);
    final LongRange r = SequencesReaderFactory.resolveRange(preReadDir, new LongRange(startId, endId));
    final boolean rename = mFlags.isSet(RENAME);
    final SequencesReader read;
    try {
      read = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(preReadDir, r);
    } catch (final FileNotFoundException e) {
      throw new NoTalkbackSlimException(ErrorType.FILE_NOT_FOUND, preReadDir.toString());
    }

    final LineWriter p;
    final File outputFile;
    if (CommonFlags.isStdio(output)) {
      outputFile = null;
      p = new LineWriter(new OutputStreamWriter(out));
    } else {
      final boolean gzip = !mFlags.isSet(NO_GZIP);
      outputFile = gzip ? new File(output + FileUtils.GZ_SUFFIX) : new File(output);
      p = new LineWriter(new OutputStreamWriter(FileUtils.createOutputStream(outputFile, gzip, false)));
    }

    try {
      final int def = mFlags.isSet(DEFAULT_QUALITY) ? (Integer) mFlags.getValue(DEFAULT_QUALITY) + (int) '!' : -1;
      process(read, p, rename, (Integer) mFlags.getValue(LINE_FLAG), def);
    } catch (InvalidParamsException e) {
      p.close();
      warnIfDeletionFails(outputFile);
      throw e;
    } finally {
      read.close();
      p.close();
    }
    return r;
  }

  static void process(final SequencesReader read, final LineWriter out, final boolean rename, final int lineLength, final int def) throws IOException, InvalidParamsException {
    if (lineLength < 0) {
      throw new IllegalArgumentException();
    }
    while (read.nextSequence()) {
      final StringBuilder name = new StringBuilder("");
      if (rename || !read.hasNames()) {
        name.append(read.currentSequenceId());
      } else {
        name.append(read.currentFullName());
      }
      out.writeln("@" + name.toString());

      final byte[] buff = new byte[(int) read.maxLength()];
      read.readCurrent(buff);
      if (lineLength == 0) {
        out.writeln(DnaUtils.bytesToSequenceIncCG(buff, 0, read.currentLength()));
        out.writeln("+" + name.toString());
        out.writeln(getScore(read, def));
      } else {
        final String dna = DnaUtils.bytesToSequenceIncCG(buff, 0, read.currentLength());
        final int dnaLen = dna.length();
        for (long i = 0; i < dnaLen; i += lineLength) {
          out.writeln(dna.substring((int) i, (int) Math.min(i + lineLength, dnaLen)));
        }
        out.writeln("+" + name.toString());
        final String qual = getScore(read, def);
        final int qualLen = qual.length();
        for (long i = 0; i < qualLen; i += lineLength) {
          out.writeln(qual.substring((int) i, (int) Math.min(i + lineLength, qualLen)));
        }
      }
    }
  }

  private static String getScore(final SequencesReader read, final int c) throws IOException, InvalidParamsException {
    if (read.hasQualityData()) {
      final byte[] quality = new byte[read.currentLength()];
      read.readCurrentQuality(quality);
      return getQuality(quality);
    } else if (c >= (int) '!') {
      if (sScore == null) {
        sScore = getSimulatedScore(read.maxLength(), (char) c);
      }
      return sScore.substring(0, read.currentLength());
    } else {
      throw new InvalidParamsException(ErrorType.INFO_ERROR, "The input SDF does not have quality data and no default was provided.");
    }
  }

  /**
   * Converts an array of bytes into a sanger-encoded quality string
   *
   * @param quality buffer containing input qualities
   * @return the quality string
   */
  public static String getQuality(byte[] quality) {
    final StringBuilder b = new StringBuilder();
    for (final byte q : quality) {
      b.append((char) (q + (byte) '!'));
    }
    return b.toString();
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  /**
   * Initialise a flags object
   * @param flags the flags object to initialise
   */
  public void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Converts SDF data into FASTQ file(s).");
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('i', INPUT, File.class, "SDF", "SDF containing sequences").setCategory(INPUT_OUTPUT);
    flags.registerRequired('o', OUTPUT_FILE, File.class, "FILE", "output filename (extension added if not present)").setCategory(INPUT_OUTPUT);
    flags.registerOptional('l', LINE_FLAG, Integer.class, "INT", "maximum number of nucleotides or amino acids to print on a line of fastq output. A value of 0 indicates no limit", 0).setCategory(UTILITY);
    flags.registerOptional('R', RENAME, "rename the reads to their consecutive number; name of first read in file is '0'").setCategory(UTILITY);
CommonFlags.initNoGzip(flags);
    flags.registerOptional('q', DEFAULT_QUALITY, Integer.class, "INT", "default quality value to use if the SDF does not contain quality data (0-63)").setCategory(UTILITY);
    flags.registerOptional(START_SEQUENCE, Long.class, "INT", "inclusive lower bound on sequence id").setCategory(CommonFlagCategories.FILTERING);
    flags.registerOptional(END_SEQUENCE, Long.class, "INT", "exclusive upper bound on sequence id").setCategory(CommonFlagCategories.FILTERING);
    flags.setValidator(VALIDATOR);
  }

  /**
   * Simulates a score
   * @param maxLength the length of the score
   * @param c the character
   * @return the score
   */
  public static String getSimulatedScore(final long maxLength, final char c) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < maxLength; i++) {
      sb.append(c);
    }
    return sb.toString();
  }
}
