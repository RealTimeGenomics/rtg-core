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

package com.rtg.sam;

import java.io.IOException;

import com.rtg.util.PopulatorFactory;
import com.rtg.util.intervals.RangeList;

import htsjdk.samtools.SAMFileHeader;

/**
 * Used to prevent having to create a separate ThreadedMultifileIterator for each reference
 * sequence being processed.
 */
public class ThreadedMultifileIteratorWrapper<T extends ReaderRecord<T> & MateInfo> implements RecordIterator<T> {

  private final RecordIterator<T> mIterator;
  private final SamReadingContext mContext;
  private T mNext;
  private int mSequenceId = -1;
  private RangeList<String> mCurrentRangeList;

  /**
   * Constructor
   * @param context the SAM reading context
   * @param populatorFactory populator factory for records
   * @throws IOException if an IO error occurs
   */
  public ThreadedMultifileIteratorWrapper(SamReadingContext context, PopulatorFactory<T> populatorFactory) throws IOException {
    final ThreadedMultifileIterator<T> it = new ThreadedMultifileIterator<>(context, populatorFactory);
    mContext = context;
    mIterator = context.filterParams().findAndRemoveDuplicates() ? new DedupifyingRecordIterator<>(it) : it;
    if (mIterator.hasNext()) {
      mNext = mIterator.next();
    } else {
      mNext = null;
    }
  }


  /**
   * Set the current reference sequence index
   * @param index reference sequence index
   */
  public void setSequenceId(int index) {
    assert index > mSequenceId : index + " > " + mSequenceId; //incoming index should always be greater than previous index
    mSequenceId = index;
    while (mNext != null && index != mNext.getSequenceId() && index > mNext.getSequenceId()) { //update mNext until we get to next sequence
      if (mIterator.hasNext()) {
        mNext = mIterator.next();
      } else {
        mNext = null;
      }
    }
    mCurrentRangeList = mContext.hasRegions() ? mContext.rangeList(mSequenceId) : null;
  }

  @Override
  public boolean hasNext() {
    if (mSequenceId < 0) {
      throw new IllegalStateException("reference is not set");
    }
    return mNext != null && mNext.getSequenceId() == mSequenceId;
  }

  @Override
  public T next() {
    final T res = mNext;
    if (mIterator.hasNext()) {
      mNext = mIterator.next();
    } else {
      mNext = null;
    }
    return res;
  }

  @Override
  public void remove() {
  }

  @Override
  public void close() throws IOException {
    mIterator.close();
  }

  @Override
  public SAMFileHeader header() {
    return mIterator.header();
  }

  @Override
  public long getTotalNucleotides() {
    return mIterator.getTotalNucleotides();
  }

  @Override
  public long getInvalidRecordsCount() {
    return mIterator.getInvalidRecordsCount();
  }

  @Override
  public long getFilteredRecordsCount() {
    return mIterator.getFilteredRecordsCount();
  }

  @Override
  public long getDuplicateRecordsCount() {
    return mIterator.getDuplicateRecordsCount();
  }

  @Override
  public long getOverCoverageRecordsCount() {
    return mIterator.getOverCoverageRecordsCount();
  }

  @Override
  public long getOutputRecordsCount() {
    return mIterator.getOutputRecordsCount();
  }

  @Override
  public long getTotalRecordsCount() {
    return mIterator.getTotalRecordsCount();
  }

  /**
   * Get the context used for input reading
   * @return the context
   */
  public SamReadingContext context() {
    return mContext;
  }

  /**
   * Get the ranges applying to the current region
   * @return the ranges.
   */
  public RangeList<String> getCurrentRangeList() {
    return mCurrentRangeList;
  }

}
