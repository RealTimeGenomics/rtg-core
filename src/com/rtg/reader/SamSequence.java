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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.sam.SamBamConstants;

import htsjdk.samtools.SAMRecord;

/**
 * Minimal wrapper to hold only necessary information about a <code>SAMRecord</code> for the
 * <code>SamBamSequenceDataSource</code>.
 */
@TestClass("com.rtg.reader.MappedSamBamSequenceDataSourceTest")
public class SamSequence {

  private static final byte READ_PAIRED_FLAG = 0x01;
  private static final byte FIRST_OF_PAIR_FLAG = 0x02;
  private static final byte READ_STRAND_FLAG = 0x04;
  private static final byte NOT_PRIMARY_ALIGNMENT_FLAG = 0x08;

  private final byte[] mReadBases;
  private final byte[] mBaseQualities;
  private final String mReadName;
  private final byte mFlags;

  private final int mProjectedSplitReadPosition;

  /**
   * Turn a <code>SAMRecord</code> into a <code>SamSequence</code>.
   * @param record the <code>SAMRecord</code> to convert.
   */
  public SamSequence(SAMRecord record) {
    assert record != null;
    //Done this way to not keep entire string of SAMRecord just for the name.
    mReadName = new String(record.getReadName().toCharArray());
    mFlags = getFlags(record);
    mReadBases = record.getReadBases();
    mBaseQualities = record.getBaseQualities();

    mProjectedSplitReadPosition = record.getAlignmentStart() * ((record.getFlags() & SamBamConstants.SAM_MATE_IS_REVERSE) != 0 ? 1 : -1);
  }

  private SamSequence(SamSequence copy, String name) {
    mReadName = name;
    mFlags = copy.mFlags;
    mReadBases = copy.mReadBases;
    mBaseQualities = copy.mBaseQualities;

    mProjectedSplitReadPosition = copy.mProjectedSplitReadPosition;
  }

  SamSequence rename(String newName) {
    return new SamSequence(this, newName);
  }

  private static byte getFlags(SAMRecord record) {
    byte flags = 0;
    if (record.getReadPairedFlag()) {
      flags ^= READ_PAIRED_FLAG;
      if (record.getFirstOfPairFlag()) {
        flags ^= FIRST_OF_PAIR_FLAG;
      }
    }
    if (record.getReadNegativeStrandFlag()) {
      flags ^= READ_STRAND_FLAG;
    }
    if (record.getNotPrimaryAlignmentFlag()) {
      flags ^= NOT_PRIMARY_ALIGNMENT_FLAG;
    }
    return flags;
  }

  /**
   * Return true if record was paired.
   * @return true if paired record.
   */
  public boolean getReadPairedFlag() {
    return (mFlags & READ_PAIRED_FLAG) != 0;
  }

  /**
   * Return true if record was first of a pair of records.
   * @return true if record was first of a pair of records.
   */
  public boolean getFirstOfPairFlag() {
    return (mFlags & FIRST_OF_PAIR_FLAG) != 0;
  }

  /**
   * Return true if record was paired.
   * @return true if paired record.
   */
  public boolean getReadNegativeStrandFlag() {
    return (mFlags & READ_STRAND_FLAG) != 0;
  }

  /**
   * Get the read bases.
   * @return the read bases.
   */
  public byte[] getReadBases() {
    return mReadBases;
  }

  /**
   * Get the base quality bytes.
   * @return the base quality bytes.
   */
  public byte[] getBaseQualities() {
    return mBaseQualities;
  }

  /**
   * Get the read name.
   * @return the read name.
   */
  public String getReadName() {
    return mReadName;
  }

  /**
   * Get read base length.
   * @return number of bases in the read.
   */
  public int getReadLength() {
      return mReadBases.length;
  }

  /**
   * @return the estimated position (0 based), or -1 if none
   */
  public int getProjectedPosition() {
    return mProjectedSplitReadPosition;
  }
}
