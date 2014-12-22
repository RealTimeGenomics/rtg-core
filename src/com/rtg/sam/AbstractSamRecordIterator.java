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
package com.rtg.sam;

import com.reeltwo.jumble.annotations.TestClass;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

/**
 * Implements the minimal header tracking and record counting likely to be used by any lowest-level <code>SAMRecord</code> RecordIterator
 */
@TestClass("com.rtg.sam.SkipInvalidRecordsIteratorTest")
public abstract class AbstractSamRecordIterator implements RecordIterator<SAMRecord> {

  protected SAMFileHeader mHeader;
  protected long mTotalNucleotides = 0;
  protected long mRecordCount = 0;
  protected long mNumInvalidRecords = 0;
  protected long mNumFilteredRecords = 0;
  protected long mNumDuplicateRecords = 0;

  /**
   * Constructor
   * @param header give it the damn header
   */
  public AbstractSamRecordIterator(SAMFileHeader header) {
    mHeader = header;
  }

  @Override
  public long getTotalRecordsCount() {
    return mRecordCount;
  }

  @Override
  public long getTotalNucleotides() {
    return mTotalNucleotides;
  }

  @Override
  public long getInvalidRecordsCount() {
    return mNumInvalidRecords;
  }

  @Override
  public long getFilteredRecordsCount() {
    return mNumFilteredRecords;
  }

  @Override
  public long getDuplicateRecordsCount() {
    return mNumDuplicateRecords;
  }

  @Override
  public long getOutputRecordsCount() {
    return getTotalRecordsCount() - getInvalidRecordsCount() - getDuplicateRecordsCount() - getFilteredRecordsCount();
  }

  @Override
  public SAMFileHeader header() {
    return mHeader;
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException("Not supported yet.");
  }
}
