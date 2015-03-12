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

import static com.rtg.util.cli.CommonFlagCategories.FILTERING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collection;
import java.util.Map;

import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.LoggedCli;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.LogStream;

/**
 * Pulls out a subset of sequences from one SDF into another SDF.
 *
 */
public final class SdfSubset extends LoggedCli {

  private static final String INPUT_FLAG = "input";
  private static final String OUTPUT_FLAG = "output";
  private static final String ID_FILE_FLAG = "id-file";
  private static final String NAMES_FLAG = "names";

  private static final Validator VALIDATOR = new Validator() {

    @Override
    public boolean isValid(CFlags flags) {
      if (!CommonFlags.validateOutputDirectory((File) flags.getValue(OUTPUT_FLAG))) {
        return false;
      }
      if (!flags.isSet(ID_FILE_FLAG) && !flags.getAnonymousFlag(0).isSet()) {
        flags.setParseMessage("Sequence ids must be specified, either explicitly, or using --" + ID_FILE_FLAG);
        return false;
      }
      return true;
    }
  };


  private byte[] mData = null;
  private byte[] mQualities = null;
  private int mWarnCount = 0;
  private long mWritten = 0;
  private SdfWriterWrapper mWriter = null;
  private SdfReaderWrapper mReader = null;
  private Map<String, Long> mNames = null;
  private final AbstractSdfWriter.SequenceNameHandler mHandler = new AbstractSdfWriter.SequenceNameHandler();

  @Override
  protected void initFlags() {
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setDescription("Extracts a subset of sequences from one SDF and outputs them to another SDF.");
    mFlags.registerRequired('i', INPUT_FLAG, File.class, "SDF", "input SDF").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, "SDF", "output SDF").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('n', NAMES_FLAG, "specify sequence names instead of sequence ids").setCategory(FILTERING);
    final Flag listFlag = mFlags.registerOptional('I', ID_FILE_FLAG, File.class, "FILE", "file containing sequence ids, or sequence names if " + NAMES_FLAG + " flag is set, one per line").setCategory(FILTERING);
    final Flag idFlag = mFlags.registerRequired(String.class, "STRING", "id of sequence to extract, or sequence name if " + NAMES_FLAG + " flag is set").setMinCount(0).setMaxCount(Integer.MAX_VALUE).setCategory(FILTERING);
    mFlags.addRequiredSet(idFlag);
    mFlags.addRequiredSet(listFlag);
    mFlags.setValidator(VALIDATOR);
  }

  @Override
  public String moduleName() {
    return "sdfsubset";
  }

  @Override
  protected File outputDirectory() {
    return (File) mFlags.getValue(OUTPUT_FLAG);
  }


  private void warnInvalidId(String seqid) {
    if (mWarnCount < 5) {
      if (mNames == null) {
        Diagnostic.warning("Invalid sequence id " + seqid + ", must be from 0 to " + (mReader.numberSequences() - 1));
      } else {
        Diagnostic.warning("Invalid sequence name " + seqid);
      }
      mWarnCount++;
      if (mWarnCount == 5) {
        Diagnostic.warning("(Only the first 5 messages shown, see log for full list.)");
      }
    } else {
      if (mNames == null) {
        Diagnostic.userLog("Invalid sequence id " + seqid);
      } else {
        Diagnostic.userLog("Invalid sequence name " + seqid);
      }
    }
  }

  private void getSequence(long seqid) throws IOException {
    if ((seqid < 0) || (seqid >= mReader.numberSequences())) {
      warnInvalidId("" + seqid);
      return;
    }
    final int length = mReader.maxLength();
    if (mData == null || mData.length < length) {
      mData = new byte[length];
      if (mReader.hasQualityData()) {
        mQualities = new byte[length];
      }
    }
    mWriter.writeSequence(mReader, seqid, mData, mQualities);
    if (++mWritten % 1000 == 0) {
      Diagnostic.progress("Extracted " + mWritten + " sequences");
    }
  }

  private void getSequences(String seqRange) throws IOException {
    if (mNames != null) {
      final Long id = mNames.get(mHandler.handleSequenceName(seqRange).label());
      if (id != null) {
        getSequence(id);
      } else {
        warnInvalidId(seqRange);
      }
      return;
    }
    final int rangePos = seqRange.indexOf("-");
    if (rangePos == -1) {
      getSequence(Long.parseLong(seqRange));
    } else {
      final String start = seqRange.substring(0, rangePos);
      final String end = seqRange.substring(rangePos + 1);
      final long startIdx;
      if (start.length() == 0) {
        startIdx = 0;
      } else {
        startIdx = Long.parseLong(start);
      }
      final long endIdx;
      if (end.length() == 0) {
        endIdx = mReader.numberSequences() - 1;
      } else {
        endIdx = Long.parseLong(end);
      }
      if (startIdx > endIdx) {
        throw new NumberFormatException("Invalid range: " + seqRange);
      }
      //System.out.println("Getting sequences from " + startIdx + " to " + endIdx);
      for (long i = startIdx; i <= endIdx; i++) {
        getSequence(i);
      }
    }
  }

  @Override
  protected int mainExec(OutputStream out, LogStream log) throws IOException {
    mReader = new SdfReaderWrapper((File) mFlags.getValue(INPUT_FLAG), false, false);

    if (mFlags.isSet(NAMES_FLAG)) {
      mNames = ReaderUtils.getSequenceNameMap(mReader.isPaired() ? mReader.left() : mReader.single());
    }
    final File output = (File) mFlags.getValue(OUTPUT_FLAG);
    mWriter = new SdfWriterWrapper(output, mReader, false);
    mWriter.setSdfId(new SdfId()); //random.nextLong());

    if (mFlags.getAnonymousFlag(0).isSet()) {
      final Collection<Object> seqs = mFlags.getAnonymousValues(0);
      for (final Object oi : seqs) {
        try {
          getSequences((String) oi);
        } catch (final NumberFormatException e) {
          warnInvalidId((String) oi);
        }
      }
    }
    if (mFlags.isSet(ID_FILE_FLAG)) {
      try (BufferedReader br = new BufferedReader(new FileReader((File) mFlags.getValue(ID_FILE_FLAG)))) {
        String line;
        while ((line = br.readLine()) != null) {
          line = line.trim();
          if (line.length() > 0) {
            try {
              getSequences(line);
            } catch (final NumberFormatException e) {
              warnInvalidId(line);
            }
          }
        }
      }
    }

    mWriter.copySourceTemplatesFile(mReader);
    mWriter.close();
    mReader.close();
    Diagnostic.progress("Extracted " + mWritten + " sequences");
    return 0;
  }

  /**
   * Main method for pulling out reads into another SDF.
   * @param args arguments for <code>SdfSubset</code>
   */
  public static void main(String[] args) {
    new SdfSubset().mainExit(args);
  }
}
