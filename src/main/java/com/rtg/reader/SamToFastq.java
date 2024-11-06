/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.reader;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.List;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommandLineFiles;
import com.rtg.launcher.CommonFlags;
import com.rtg.sam.SamFilter;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.BaseFile;
import com.rtg.util.io.FileUtils;

/**
 * Convert SAM to FASTQ files. Only works on paired end SAM/BAM.
 * Memory use can be reduced by using pre-collated (see samtools collate) files.
 */
public class SamToFastq extends AbstractCli {
  private static final int MAX_WARNINGS = 5;

  private static final String MODULE_NAME = "sam2fastq";
  private static final String INPUT_FLAG = "input";
  private static final String COLLATED = "collated";

  @Override
  public String description() {
    return "Extract fastq from paired-end SAM/BAM files";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription(description());
    initFlags(mFlags);
  }

  private static void initFlags(CFlags flags) {
    CommonFlagCategories.setCategories(flags);
    final Flag<File> input = flags.registerOptional('i', INPUT_FLAG, File.class, CommonFlags.FILE, "input file to convert to fastq").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    final Flag<File> inputList = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "list of input file to convert to fastq").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, CommonFlags.FILE, "output file base name").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional(COLLATED, "if set, assume input is collated to reduce memory use").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional('Z', CommonFlags.NO_GZIP, "Do not gzip output").setCategory(CommonFlagCategories.UTILITY);
    flags.addRequiredSet(input);
    flags.addRequiredSet(inputList);
    flags.setValidator(f -> {
      if (!f.checkXor(INPUT_FLAG, CommonFlags.INPUT_LIST_FLAG)) {
        f.setParseMessage("Must set one of --" + INPUT_FLAG + " or --" + CommonFlags.INPUT_LIST_FLAG);
        return false;
      }
      return true;
    });
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final List<File> files = new CommandLineFiles(CommonFlags.INPUT_LIST_FLAG, INPUT_FLAG).getFileList(mFlags);

    final SequenceDataSource ds = mFlags.isSet(COLLATED)
      ? new CollatedDataSource(new FileStreamIterator(files), true, false, null)
      : MappedSamBamSequenceDataSource.fromInputFiles(files, true, false, true, null);

    final BaseFile baseFile = FastqUtils.baseFile((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG), !mFlags.isSet(CommonFlags.NO_GZIP));
    boolean left = true;
    try (FastqWriter outLeft = new FastqWriter(new OutputStreamWriter(FileUtils.createOutputStream(baseFile, "_1")))) {
      try (FastqWriter outRight = new FastqWriter(new OutputStreamWriter(FileUtils.createOutputStream(baseFile, "_2")))) {
        while (ds.nextSequence()) {
          final FastqWriter w = left ? outLeft : outRight;
          w.write(ds.name(), ds.sequenceData(), ds.qualityData(), ds.currentLength());
          left = !left;
        }
        
        try (OutputStreamWriter writer = new OutputStreamWriter(out)) {
          writer.write(ds.toString());
        }
      }
    }
    return 0;
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  static class CollatedDataSource extends SamBamSequenceDataSource {

    private String mLastName;
    private FragSide mLastFragSide;
    private boolean mFragDone;

    //stats
    private long mValidRecordPairs;
    private long mValidRecords;
    private long mSkippedRecords;

    private long mBadPairs;

    private enum FragSide {
      FIRST_OF_FRAGMENT,
      LAST_OF_FRAGMENT
    }


    CollatedDataSource(FileStreamIterator inputs, boolean paired, boolean flattenPaired, SamFilter filter) {
      super(inputs, paired, flattenPaired, filter);
    }

    @Override
    protected boolean nextRecords() throws IOException {
      boolean done = false;
      while (!done) {
        final SamSequence samSequence = nextRecord();
        if (samSequence != null) {
          if (!samSequence.getReadName().equals(mLastName)) {
            if (mLastName != null && !mFragDone) {
              ++mBadPairs;
              if (mBadPairs < MAX_WARNINGS) {
                Diagnostic.warning("No pair found for sequence: " + mLastName);
              } else if (mBadPairs == MAX_WARNINGS) {
                Diagnostic.warning("No pair found for sequence: " + mLastName + " further messages of this type will be suppressed");
              }
            }
            mLastName = samSequence.getReadName();
            if (samSequence.getFirstOfPairFlag()) {
              mLastFragSide = FragSide.FIRST_OF_FRAGMENT;
            } else {
              mLastFragSide = FragSide.LAST_OF_FRAGMENT;
            }
            placePairedRecord(samSequence);
            ++mValidRecords;
            mFragDone = false;
          } else if (!mFragDone && samSequence.getFirstOfPairFlag() != (mLastFragSide == FragSide.FIRST_OF_FRAGMENT)) {
              placePairedRecord(samSequence);
              ++mValidRecords;
              ++mValidRecordPairs;
              mFragDone = true;
              done = true;
          } else {
            ++mSkippedRecords;
          }
        } else {
          break;
        }
      }
      return haveNextRecords();
    }

    @Override
    protected void checkSortOrder() {
    }

    public String toString() {
      final StringBuilder sb = new StringBuilder();
      sb.append("Record Pairs written: ").append(mValidRecordPairs);
      sb.append(StringUtils.LS);
      sb.append("Records written: ").append(mValidRecords);
      sb.append(StringUtils.LS);
      sb.append("Duplicate records skipped: ").append(mSkippedRecords);
      sb.append(StringUtils.LS);
      if (mBadPairs > 0) {
        sb.append(StringUtils.LS);
        sb.append("Arms with missing pairs: ").append(mBadPairs);
        sb.append(StringUtils.LS);
      }
      return sb.toString();
    }
  }

}
