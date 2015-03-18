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

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.UTILITY;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.mode.SequenceType;
import com.rtg.util.PortableRandom;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.InformationType;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;
import com.rtg.util.io.LogStream;

/**
 * Splits an SDF file into parts.
 *
 */
public final class SdfSplitter extends LoggedCli {

  private static final String MODULE_NAME = "sdfsplit";

  private static final String OUTPUT_FLAG = "output";
  private static final String COUNT_FLAG = "num-sequences";
  private static final String ENABLE_MEMORY_READER = "in-memory";
  private static final String DISABLE_DUPLICATE_DETECTOR = "allow-duplicate-names";
  private static final String XFORCE_COMPRESS = "Xforce-compress";

  /**
   * @return current name of the module
   */
  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  /**
   * construct a flags object
   */
  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Splits SDF data into multiple equal segments.");
    CommonFlagCategories.setCategories(mFlags);
    final Flag inFlag = mFlags.registerRequired(File.class, "SDF", "input SDF");
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    inFlag.setCategory(INPUT_OUTPUT);
    final Flag listFlag = mFlags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of SDFs (1 per line)").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, "DIR", "output base directory (must be empty if present)").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('n', COUNT_FLAG, Integer.class, "INT", "number of reads per output").setCategory(UTILITY);
    mFlags.registerOptional(DISABLE_DUPLICATE_DETECTOR, "disable checking for duplicate sequence names").setCategory(UTILITY);

    mFlags.registerOptional(ENABLE_MEMORY_READER, "process in memory (faster but requires more RAM)").setCategory(UTILITY);
    mFlags.registerOptional(XFORCE_COMPRESS, "force compression for output SDFs").setCategory(UTILITY);
    mFlags.addRequiredSet(inFlag);
    mFlags.addRequiredSet(listFlag);
    mFlags.setValidator(VALIDATOR);
  }

  private static final Validator VALIDATOR = new Validator() {
    @Override
    public boolean isValid(final CFlags flags) {

      if (!CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE)) {
        return false;
      }

      if (!CommonFlags.validateOutputDirectory((File) flags.getValue(OUTPUT_FLAG))) {
        return false;
      }

      if ((Integer) flags.getValue(COUNT_FLAG) < 1) {
        flags.setParseMessage("Value for --" + COUNT_FLAG + " must be greater than 0.");
        return false;
      }

      return true;
    }
  };

  static void split(final List<File> inDirs, final File outDir, final long sequencesPerOutput, boolean useMem, boolean disableDupeDetect, boolean forceCompression) throws IOException {
    final PortableRandom random = new PortableRandom();
    final boolean isPaired;
    final boolean hasQuality;
    final boolean hasNames;
    long guid = 0;
    final SequenceType type;
    int maxLength;
    SdfReaderWrapper reader;
    long totalSequences;
    try {
      reader = new SdfReaderWrapper(inDirs.get(0), false, true);
    } catch (final FileNotFoundException e) {
      throw new NoTalkbackSlimException(ErrorType.FILE_NOT_FOUND, inDirs.get(0).toString());
    }
    isPaired = reader.isPaired();
    hasQuality = reader.hasQualityData();
    hasNames = reader.hasNames();
    while (guid == 0) {
      guid = random.nextLong();
    }
    maxLength = reader.maxLength();
    type = reader.type();
    totalSequences = reader.numberSequences();
    reader.close();
    for (int i = 1; i < inDirs.size(); i++) {
      final File dir = inDirs.get(i);
      try {
        reader = new SdfReaderWrapper(dir, false, isPaired);
      } catch (final FileNotFoundException e) {
        throw new NoTalkbackSlimException(ErrorType.FILE_NOT_FOUND, dir.toString());
      }
      if (reader.isPaired() != isPaired) {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Mixed paired and non-paired SDF files provided.");
      }
      if (reader.hasQualityData() != hasQuality) {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Mixed quality and non-quality SDF files provided.");
      }
      if (reader.hasNames() != hasNames) {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Mixed names and non-names SDF files provided.");
      }
      if (!reader.type().equals(type)) {
        throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Mixed DNA and Protein SDF files provided.");
      }
      maxLength = Math.max(maxLength, reader.maxLength());
      totalSequences += reader.numberSequences();
      reader.close();
    }

    long numSequences = 0;
    long numOutputs = 0;
    final byte[] data = new byte[maxLength];
    final byte[] quality = new byte[maxLength];
    SdfWriterWrapper writer = null;
    final NameDuplicateDetector dupDetector = disableDupeDetect || !hasNames ? NameDuplicateDetector.getNullDetector() : new NameDuplicateDetector(totalSequences);
    reader = null;
    int readerNumber = -1;
    final int numberOutputs = (int) Math.ceil(totalSequences / (double) sequencesPerOutput);
    Diagnostic.info(InformationType.INFO_USER, "sdfsplit is creating " + numberOutputs + " SDF" + (numberOutputs != 1 ? "s" : ""));
    try {
      for (final File dir : inDirs) {
        try {
          reader = new SdfReaderWrapper(dir, useMem, false);
          readerNumber++;
        } catch (final FileNotFoundException e) {
          throw new NoTalkbackSlimException(ErrorType.FILE_NOT_FOUND, dir.toString());
        }
        for (long seq = 0; seq < reader.numberSequences(); seq++) {
          if (writer == null) {
            final String fname;
            fname = String.format("%06d", numOutputs);
            writer = new SdfWriterWrapper(new File(outDir, fname), reader, forceCompression);
            writer.setSdfId(new SdfId()); //guid + numOutputs);
            numOutputs++;
            if (inDirs.size() == 1) {
              writer.copySourceTemplatesFile(reader);
            }
          }
          dupDetector.addPair(reader.name(seq), (int) seq, readerNumber);
          writer.writeSequence(reader, seq, data, quality);
          if (++numSequences == sequencesPerOutput) {
            writer.close();
            writer = null;
            numSequences = 0;
          }
        }
        reader.close();
        reader = null;
      }
      if (!disableDupeDetect) {
        readerNumber = 0;
        final SequencesReader[] readers = new SequencesReader[inDirs.size()];
        final SdfReaderWrapper[] readerWrappers = new SdfReaderWrapper[inDirs.size()];
        try {
          for (final File dir : inDirs) {
            try {
              readerWrappers[readerNumber] = new SdfReaderWrapper(dir, false, false);
            } catch (final FileNotFoundException e) {
              throw new NoTalkbackSlimException(ErrorType.FILE_NOT_FOUND, dir.toString());
            }
            readers[readerNumber] = readerWrappers[readerNumber].isPaired() ? readerWrappers[readerNumber].left() : readerWrappers[readerNumber].single();
            readerNumber++;
          }
          if (dupDetector.checkSequenceDuplicates(readers, new File(outDir.getPath(), "duplicate-names.txt"))) {
            Diagnostic.warning(WarningType.INFO_WARNING, "Duplicate Sequence Names in Input");
            //throw new NoTalkbackSlimException(ErrorType.INFO_ERROR, "Duplicate Sequence Names in Input");
          }
        } finally {
          for (final SdfReaderWrapper r : readerWrappers) {
            if (r != null) {
              try {
                r.close();
              } finally { }
            }
          }
        }
      }
    } finally {
      try {
        if (reader != null) {
          reader.close();
        }
      } finally {
        if (writer != null) {
          writer.close();
        }
      }
    }
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    final File outDir = outputDirectory();
    final int count = (Integer) mFlags.getValue(COUNT_FLAG);
    final boolean useMem = mFlags.isSet(ENABLE_MEMORY_READER);
    final List<File> inDirs = CommonFlags.getFileList(mFlags, CommonFlags.INPUT_LIST_FLAG, null, true);
    final boolean forceCompression = mFlags.isSet(XFORCE_COMPRESS);
    split(inDirs, outDir, count, useMem, mFlags.isSet(DISABLE_DUPLICATE_DETECTOR), forceCompression);
    return 0;
  }

  /**
   * Splits an SDF file into parts.
   *
   * @param args Command line arguments
   */
  public static void main(final String[] args) {
    new SdfSplitter().mainExit(args);
  }
}

