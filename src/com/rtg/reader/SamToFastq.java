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
 * Only works on paired end pre-collated (see samtools collate) SAM/BAM
 */
public class SamToFastq extends AbstractCli {
  private static final int MAX_WARNINGS = 5;

  private static final String MODULE_NAME = "sam2fastq";
  private static final String INPUT_FLAG = "input";



  @Override
  public String description() {
    return "Extract fastq from pre-collated paired-end SAM/BAM files";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription(description());
    initFlags(mFlags);
  }

  private static void initFlags(CFlags flags) {
    CommonFlagCategories.setCategories(flags);
    final Flag<File> input = flags.registerOptional('i', INPUT_FLAG, File.class, CommonFlags.FILE, "input file to convert to fastq").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    final Flag<File> inputList = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, CommonFlags.FILE, "List of input file tos convert to fastq").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, CommonFlags.FILE, "output file base name").setCategory(CommonFlagCategories.INPUT_OUTPUT);
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
    final CollatedDataSource ds = new CollatedDataSource(new FileStreamIterator(files), true, false, null);

    final BaseFile baseOutFile = FastqUtils.baseFile((File) mFlags.getValue(CommonFlags.OUTPUT_FLAG), !mFlags.isSet(CommonFlags.NO_GZIP));
    final File outLeftName = baseOutFile.suffixedFile("_1");
    final File outRightName = baseOutFile.suffixedFile("_2");
    boolean left = true;
    try (FastqWriter outLeft = new FastqWriter(new OutputStreamWriter(FileUtils.createOutputStream(outLeftName)))) {
      try (FastqWriter outRight = new FastqWriter(new OutputStreamWriter(FileUtils.createOutputStream(outRightName)))) {
        while (ds.nextSequence()) {
          final FastqWriter w = left ? outLeft : outRight;
          w.write(ds.name(), ds.sequenceData(), ds.qualityData(), ds.currentLength());
          left = !left;
        }
        try (OutputStreamWriter writer = new OutputStreamWriter(out)) {
          writer.write("Record Pairs written: " + ds.mValidRecordPairs);
          writer.write(StringUtils.LS);
          writer.write("Records written: " + ds.mValidRecords);
          writer.write(StringUtils.LS);
          writer.write("Duplicate records skipped: " + ds.mSkippedRecords);
          writer.write(StringUtils.LS);
          if (ds.mBadPairs > 0) {
            writer.write(StringUtils.LS);
            writer.write("Arms with missing pairs: " + ds.mBadPairs);
            writer.write(StringUtils.LS);
          }
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
  }

}
