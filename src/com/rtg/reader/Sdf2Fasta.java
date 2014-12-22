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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.Locale;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.ProteinFastaSymbolTable;
import com.rtg.mode.SequenceType;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.io.FileUtils;

/**
 * Convert a preread dataset into FASTA format.
 *
 */
public class Sdf2Fasta extends AbstractCli {

  static final String MODULE_NAME = "sdf2fasta";

  static final String FINPUT_FLAG = "input";
  static final String FOUTPUT_FLAG = "output";
  static final String LINE_FLAG = "line-length";
  static final String RENAME = Sdf2Fastq.RENAME;
  static final String CHUNKS_FLAG = "Xchunks";
  static final String COLORSPACE_FLAG = "Xcolorspace";

  static final String START_SEQUENCE = "start-id";
  static final String END_SEQUENCE = "end-id";

  static final String FASTA_EXT_LGE = ".fasta";
  static final String FASTA_EXT_SRT = ".fa";

  /**
   * @return current name of the module
   */
  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  static BufferedWriter getStream(final String name, final boolean append, boolean gzip) throws IOException {
    if (CommonFlags.isStdio(name)) {
      return new BufferedWriter(new OutputStreamWriter(FileUtils.getStdoutAsOutputStream()), FileUtils.BUFFERED_STREAM_SIZE);
    }
    return new BufferedWriter(new OutputStreamWriter(FileUtils.createOutputStream(gzip ? new File(name + FileUtils.GZ_SUFFIX) : new File(name), gzip, append)), FileUtils.BUFFERED_STREAM_SIZE);
  }

  private static final byte[] DNA_MAPPING = new DNAFastaSymbolTable().getOrdinalToAsciiTable();
  private static final byte[] PROTEIN_MAPPING = new ProteinFastaSymbolTable().getOrdinalToAsciiTable();
  private static final byte[] COLORSPACE = "0123103223013210".getBytes();

  static boolean toColorspace(final byte[] residues, final int valid) {
    if (valid > 0) {
      int current = residues[0] - 1;
      if (current == -1) {
        // Unknown nucleotide at start of sequence
        return false;
      }
      residues[0] = DNA_MAPPING[residues[0]];
      for (int k = 1; k < valid; k++) {
        final int t = residues[k] - 1;
        if (t == -1) {
          // Unknown nucleotide detected
          return false;
        }
        final byte v = COLORSPACE[(current << 2) + t];
        current = t;
        residues[k] = v;
      }
    }
    return true;
  }

  private int mChunks = 1; // Maximum number of output files
  private int mLineLength = 0; // Maximum residues per line -- 0 denotes infinite line length
  private boolean mColorspace = false;  // Output in colorspace format
  private boolean mGzip = false; // Compress output
  private boolean mRename = false; // Replace sequence names with their SDF sequence ID

  /**
   * Convert preread information into FASTA format.  Will write files starting
   * with the given prefix and ending with <code>.fasta</code>.
   *
   * @param input preread source
   * @param prefix filename root
   * @param postfix for paired end
   * @throws IOException if an I/O error occurs
   * @throws IllegalArgumentException if <code>chunks</code> is less than 1 or <code>lineLength</code> is negative.
   * @throws NullPointerException if <code>input</code> or <code>prefix</code> is null.
   */
  public void doPrereadToFasta(final SequencesReader input, final String prefix, final String postfix) throws IOException {
    if (mChunks < 1 || mLineLength < 0) {
      throw new IllegalArgumentException();
    }
    if (prefix == null) {
      throw new NullPointerException();
    }
    final String prefixPath = new File(prefix).getAbsolutePath();
    final long maxLength = input.maxLength();
    if (maxLength > Integer.MAX_VALUE) {
      throw new NoTalkbackSlimException(ErrorType.SEQUENCE_LENGTH_ERROR);
    }
    final byte[] residues = new byte[(int) maxLength];
    final byte[] mapping;
    if (input.type() == SequenceType.DNA) {
      mapping = DNA_MAPPING;
    } else {
      if (mColorspace) {
        throw new NoTalkbackSlimException(ErrorType.COLORSPACE_DNA);
      }
      mapping = PROTEIN_MAPPING;
    }

    try {
      final BufferedWriter[] outputs = new BufferedWriter[mChunks];
      outputs[0] = mChunks == 1 ? getStream(CommonFlags.isStdio(prefix) ? STDIO_NAME : (prefixPath + postfix), false, mGzip) : null;
      try {
        long containedUnknowns = 0;
        //final long totalLength = input.totalLength();
        //final PrereadToFasta p2f = new PrereadToFasta(); // for progress
        int c = 0; // file number to write to (if relevant)
        while (input.nextSequence()) {
          if (outputs[c] == null) {
            // i.e. have multiple mChunks
            outputs[c] = getStream(CommonFlags.isStdio(prefix) ? STDIO_NAME : (prefixPath + c + postfix), true, mGzip);
          }
          final int valid = input.readCurrent(residues);
          // replace internal codes by ASCII bytes in situ
          if (mColorspace) {
            if (!toColorspace(residues, valid)) {
              containedUnknowns++;
              continue;
            }
          } else {
            for (int i = 0; i < valid; i++) {
              residues[i] = mapping[residues[i]];
            }
          }
          outputs[c].write(">" + (!mRename && input.hasNames() ? input.currentFullName() : ("" + input.currentSequenceId())));
          outputs[c].newLine();
          // write actual ASCII bytes to output
          if (mLineLength == 0) {
            outputs[c].write(new String(residues, 0, valid));
            outputs[c].newLine();
          } else {
            for (long k = 0; k < valid; k += mLineLength) {
              outputs[c].write(new String(residues, (int) k, Math.min(mLineLength, valid - (int) k)));
              outputs[c].newLine();
            }
          }
          c = ++c % mChunks;
        }
        if (containedUnknowns > 0) {
          Diagnostic.warning(WarningType.COLORSPACE_WARNING, String.valueOf(containedUnknowns));
        }
      } finally {
        IOException ex = null;
        for (final BufferedWriter out : outputs) {
          try {
            if (out != null) {
              out.close();
            }
          } catch (final IOException e) {
            if (ex == null) {
              ex = e;
            }
          }
        }
        if (ex != null) {
          throw ex;
        }
      }
    } catch (final IOException e) {
      // Ignore broken pipe error so we don't die on PrereadToFasta | head etc.
      if (!e.getMessage().contains("Broken pipe")) {
        throw e;
      }
    }
  }



  private static final Validator VALIDATOR = new Validator() {
    @Override
    public boolean isValid(final CFlags flags) {
      if ((Integer) flags.getValue(CHUNKS_FLAG) < 1) {
        Diagnostic.error(ErrorType.EXPECTED_POSITIVE, CHUNKS_FLAG);
        return false;
      } else if ((Integer) flags.getValue(CHUNKS_FLAG) > 400) {
        Diagnostic.error(ErrorType.BAD_FILE_CHUNK_SIZE);
        return false;
      }
      if ((Integer) flags.getValue(LINE_FLAG) < 0) {
        Diagnostic.error(ErrorType.EXPECTED_NONNEGATIVE, LINE_FLAG);
        return false;
      }
      final File input = (File) flags.getValue(FINPUT_FLAG);
      if (!CommonFlags.validateSDF(input)) {
        return false;
      }
      return CommonFlags.validateStartEnd(flags, START_SEQUENCE, END_SEQUENCE);
    }
  };

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  /**
   * initialise a flags object
   * @param flags the flags object to initialise
   */
  public void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Converts SDF data into FASTA file(s).");
    CommonFlagCategories.setCategories(flags);
    flags.registerRequired('i', FINPUT_FLAG, File.class, "SDF", "SDF containing sequences").setCategory(INPUT_OUTPUT);
    flags.registerRequired('o', FOUTPUT_FLAG, String.class, "FILE", "output filename (extension added if not present)").setCategory(INPUT_OUTPUT);

    flags.registerOptional('c', CHUNKS_FLAG, Integer.class, "INT", "number of output files to write. The output will spread into at most this many files. Should be positive (400 maximum)", 1).setCategory(INPUT_OUTPUT);
    flags.registerOptional('l', LINE_FLAG, Integer.class, "INT", "maximum number of nucleotides or amino acids to print on a line of fasta output. A value of 0 indicates no limit", 0).setCategory(UTILITY);
CommonFlags.initNoGzip(flags);
    flags.registerOptional('R', RENAME, "rename the reads to their consecutive number; name of first read in file is '0'").setCategory(UTILITY);
    flags.registerOptional(COLORSPACE_FLAG, "write SOLiD colorspace format results. This option will only work on short sequences. Sequences containing unknown residues are skipped and their count reported at the end of the run. This option can only be used with DNA").setCategory(UTILITY);
    flags.registerOptional(START_SEQUENCE, Long.class, "INT", "inclusive lower bound on sequence id").setCategory(CommonFlagCategories.FILTERING);
    flags.registerOptional(END_SEQUENCE, Long.class, "INT", "exclusive upper bound on sequence id").setCategory(CommonFlagCategories.FILTERING);
    flags.setValidator(VALIDATOR);
  }

  /**
   * Convert an SDF to FASTA
   *
   * @param args command line arguments
   */
  public static void main(final String[] args) {
    new Sdf2Fasta().mainExit(args);
  }

  @Override
  protected int mainExec(final OutputStream out, final PrintStream err) throws IOException {
    final PrintStream outStream = new PrintStream(out);
    try {
      mChunks = (Integer) mFlags.getValue(CHUNKS_FLAG);
      mLineLength = (Integer) mFlags.getValue(LINE_FLAG);
      mColorspace = mFlags.isSet(COLORSPACE_FLAG);
      mGzip = !mFlags.isSet(NO_GZIP);
      mRename = mFlags.isSet(RENAME);

      final File input = (File) mFlags.getValue(FINPUT_FLAG);

      // Strip any FASTA like suffix from the file name
      String output = (String) mFlags.getValue(FOUTPUT_FLAG);
      String ext = FASTA_EXT_LGE;
      if (FileUtils.isGzipFilename(output)) {
        output = output.substring(0, output.length() - FileUtils.GZ_SUFFIX.length());
      }
      if (output.toLowerCase(Locale.getDefault()).endsWith(FASTA_EXT_LGE)) {
        ext = output.substring(output.length() - FASTA_EXT_LGE.length());
        output = output.substring(0, output.length() - FASTA_EXT_LGE.length());
      } else if (output.toLowerCase(Locale.getDefault()).endsWith(FASTA_EXT_SRT)) {
        ext = output.substring(output.length() - FASTA_EXT_SRT.length());
        output = output.substring(0, output.length() - FASTA_EXT_SRT.length());
      }
      final boolean isPaired = ReaderUtils.isPairedEndDirectory(input);
      if (isPaired) {
        if (CommonFlags.isStdio(output)) {
          err.println("Sending paired-end data to stdout is not supported.");
          return 1;
        }
        final LongRange calculatedRegion = doPrereadToFasta(ReaderUtils.getLeftEnd(input), output, "_1" + ext , null);
        doPrereadToFasta(ReaderUtils.getRightEnd(input), output, "_2" + ext, calculatedRegion);
      } else {
        doPrereadToFasta(input, output, ext, null);
      }
      Diagnostic.deleteLog(); // was successful execution
      return 0;
    } finally { //make sure listener removed so later unit tests not compromised
      outStream.flush();
    }
  }

  //Calculated region for sloppy end
  private LongRange doPrereadToFasta(File input, String prefix, String postfix, LongRange calculatedRegion) throws IOException {
    final long startId = mFlags.isSet(START_SEQUENCE) ? (Long) mFlags.getValue(START_SEQUENCE) : LongRange.MISSING;
    final long endId = calculatedRegion != null ? calculatedRegion.getEnd() : (mFlags.isSet(END_SEQUENCE) ? (Long) mFlags.getValue(END_SEQUENCE) : LongRange.MISSING);
    final LongRange r = SequencesReaderFactory.resolveRange(input, new LongRange(startId, endId));
    SequencesReader reader;
    try {
      reader = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(input, r);
    } catch (final FileNotFoundException e) {
      // Slightly better I/O reporting than the default provided by AbstractCli
      if (input.isDirectory()) {
        throw new NoTalkbackSlimException(ErrorType.SDF_INDEX_NOT_VALID, input.toString());
      } else if (input.exists()) {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "The specified file, \"" + input.getPath() + "\", is not an SDF.");
      } else {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "The specified SDF, \"" + input.getPath() + "\", does not exist.");
      }
    }
    try {
      doPrereadToFasta(reader, prefix, postfix);
    } finally {
      reader.close();
    }
    return r;
  }
}

