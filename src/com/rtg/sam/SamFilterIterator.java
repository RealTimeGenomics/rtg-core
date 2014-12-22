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

import java.io.IOException;
import java.util.NoSuchElementException;

import com.reeltwo.jumble.annotations.TestClass;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

/**
 * Sam iterator that supports filtering
 */
@TestClass(value = {"com.rtg.variant.cnv.CnvProductTaskTest"})
public class SamFilterIterator implements RecordIterator<SAMRecord> {

  private final RecordIterator<SAMRecord> mBaseIterator;

  private final SamFilter mSamFilter;

  private SAMRecord mNextRecord;

  private long mFilteredRecords;
  private long mValidCount = 0;

  private boolean mIsClosed;

  /**
   * Constructor
   * @param iterator the iterator to filter records from.
   * @param filter the sam record filter
   */
  public SamFilterIterator(RecordIterator<SAMRecord> iterator, SamFilter filter) {
    mBaseIterator = iterator;
    mSamFilter = filter;
    mFilteredRecords = 0;
  }

  @Override
  public void close() throws IOException {
    if (!mIsClosed) {
      mBaseIterator.close();
      mIsClosed = true;
    }
  }

  @Override
  public long getTotalNucleotides() {
    return mBaseIterator.getTotalNucleotides();
  }

  @Override
  public long getInvalidRecordsCount() {
    return mBaseIterator.getInvalidRecordsCount();
  }

  @Override
  public long getDuplicateRecordsCount() {
    return mBaseIterator.getDuplicateRecordsCount();
  }

  @Override
  public long getFilteredRecordsCount() {
    return mFilteredRecords + mBaseIterator.getFilteredRecordsCount();
  }

  @Override
  public long getOutputRecordsCount() {
    return mValidCount;
  }

  @Override
  public long getTotalRecordsCount() {
    return mBaseIterator.getTotalRecordsCount();
  }

  private void pullNextRecord() {
    if (!mBaseIterator.hasNext()) {
      mNextRecord = null;
      return;
    }
    mNextRecord = mBaseIterator.next();
  }

  @Override
  public boolean hasNext() {
    if (mNextRecord == null) {
      while (true) {
        pullNextRecord();
        if (mNextRecord != null) {
          if (!mSamFilter.acceptRecord(mNextRecord)) {
            mFilteredRecords++;
            continue;
          }
          break;
        } else {
          break;
        }
      }
    }
    return mNextRecord != null;
  }

  @Override
  public SAMRecord next() {
    if (!hasNext()) {
      throw new NoSuchElementException();
    }
    final SAMRecord ret = mNextRecord;
    mNextRecord = null;
    mValidCount++;
    return ret;
  }

  @Override
  public void remove() {
    throw new UnsupportedOperationException("Not supported.");
  }

  @Override
  public SAMFileHeader header() {
    return mBaseIterator.header();
  }

}
