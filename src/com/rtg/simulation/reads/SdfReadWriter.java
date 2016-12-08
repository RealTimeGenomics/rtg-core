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
package com.rtg.simulation.reads;

import java.io.File;
import java.io.IOException;

import com.rtg.mode.SequenceType;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SdfId;
import com.rtg.reader.SdfWriter;
import com.rtg.reader.SourceTemplateReadWriter;
import com.rtg.util.Constants;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.io.FileUtils;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 * Write simulated reads to SDF.
 */
public class SdfReadWriter implements ReadWriter {

  private final boolean mIsPaired;
  private final SdfWriter mLeft;
  private final SdfWriter mRight;
  private final SdfWriter mSingle;
  private int mTotal = 0;
  private SdfId[] mTemplateSetIds = null;
  private SdfId mOriginalReference = null;
  private boolean mExpectLeft = true;

  /**
   * Constructor
   * @param outputDir directory to create SDF in.
   * @param isPaired whether SDF should be paired
   * @param type SDF type
   * @param names true if we should write read names
   * @param quality true if we should write quality
   */
  public SdfReadWriter(File outputDir, boolean isPaired, PrereadType type, boolean names, boolean quality) {
    mIsPaired = isPaired;
    if (mIsPaired) {
      FileUtils.ensureOutputDirectory(outputDir);
      mSingle = null;
      mLeft = new SdfWriter(new File(outputDir, "left"), Constants.MAX_FILE_SIZE, type, quality, names, true, SequenceType.DNA);
      mLeft.setPrereadArm(PrereadArm.LEFT);
      mLeft.setCommandLine(CommandLine.getCommandLine());
      mRight = new SdfWriter(new File(outputDir, "right"), Constants.MAX_FILE_SIZE, type, quality, names, true, SequenceType.DNA);
      mRight.setPrereadArm(PrereadArm.RIGHT);
      mRight.setCommandLine(CommandLine.getCommandLine());
      //mRight.setOldSdfId(mLeft.getOldSdfId());
      mRight.setSdfId(mLeft.getSdfId());

    } else {
      mLeft = null;
      mRight = null;
      mSingle = new SdfWriter(outputDir, Constants.MAX_FILE_SIZE, type, quality, names, true, SequenceType.DNA);
      mSingle.setPrereadArm(PrereadArm.UNKNOWN);
      mSingle.setCommandLine(CommandLine.getCommandLine());
    }
  }

  @Override
  public void identifyTemplateSet(SdfId... templateIds) {
    mTemplateSetIds = templateIds;
  }

  @Override
  public void identifyOriginalReference(SdfId referenceId) {
    mOriginalReference = referenceId;
  }

  /**
   * Set the SDF comment
   * @param comment the comment text
   */
  public void setComment(String comment) {
    if (mIsPaired) {
      mLeft.setComment(comment);
      mRight.setComment(comment);
    } else {
      mSingle.setComment(comment);
    }
  }

  /**
   * Set the SDF read group
   * @param readGroup the read group information
   */
  public void setReadGroup(SAMReadGroupRecord readGroup) {
    if (readGroup != null) {
      if (mIsPaired) {
        mLeft.setReadGroup(readGroup.toString());
        mRight.setReadGroup(readGroup.toString());
      } else {
        mSingle.setReadGroup(readGroup.toString());
      }
    }
  }

  @Override
  public void writeRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    if (mIsPaired) {
      throw new IllegalStateException();
    }
    writeSequence(mSingle, mTotal + " " + name, data, qual, length);
    ++mTotal;
  }

  @Override
  public void writeLeftRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    if (!mIsPaired) {
      throw new IllegalStateException();
    }
    if (!mExpectLeft) {
      throw new IllegalStateException();
    }
    writeSequence(mLeft, mTotal + " " + name, data, qual, length);
    mExpectLeft = !mExpectLeft;
  }

  @Override
  public void writeRightRead(String name, byte[] data, byte[] qual, int length) throws IOException {
    if (!mIsPaired) {
      throw new IllegalStateException();
    }
    if (mExpectLeft) {
      throw new IllegalStateException();
    }
    writeSequence(mRight, mTotal + " " + name, data, qual, length);
    mExpectLeft = !mExpectLeft;
    ++mTotal;
  }

  private void writeSequence(SdfWriter writer, String name, byte[] data, byte[] qual, int length) throws IOException {
    writer.startSequence(name);
    writer.write(data, qual, length);
    writer.endSequence();
  }


  /**
   * Close method for the writers.
   * @throws IOException whenever
   */
  @Override
  public void close() throws IOException {
    if (mIsPaired) {
      if (!mExpectLeft) {
        throw new IOException("Left and Right arms were not balanced during simulation!");
      }
      if (mLeft != null) {
        mLeft.close();
        SourceTemplateReadWriter.writeTemplateMappingFile(mLeft.directory(), mTemplateSetIds);
        SourceTemplateReadWriter.writeMutationMappingFile(mLeft.directory(), mOriginalReference);
      }
      if (mRight != null) {
        mRight.close();
        SourceTemplateReadWriter.writeTemplateMappingFile(mRight.directory(), mTemplateSetIds);
        SourceTemplateReadWriter.writeMutationMappingFile(mRight.directory(), mOriginalReference);
      }
    } else {
      if (mSingle != null) {
        mSingle.close();
        SourceTemplateReadWriter.writeTemplateMappingFile(mSingle.directory(), mTemplateSetIds);
        SourceTemplateReadWriter.writeMutationMappingFile(mSingle.directory(), mOriginalReference);
      }
    }
  }

  @Override
  public int readsWritten() {
    return mTotal;
  }

}
