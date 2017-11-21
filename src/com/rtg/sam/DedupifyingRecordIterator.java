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

import htsjdk.samtools.SAMFileHeader;

/**
 * Used when one wants the functionality of the duplicate removal around a {@link RecordIterator}.
 * @param <T> type for iterator
 */
public class DedupifyingRecordIterator<T extends ReaderRecord<T> & MateInfo> extends DedupifyingIterator<T> implements RecordIterator<T> {

  private final RecordIterator<T> mRecordWrapped;

  /**
   * @param it iterator to wrap
   */
  public DedupifyingRecordIterator(RecordIterator<T> it) {
    super(it);
    mRecordWrapped = it;
  }

  @Override
  public long getDuplicateRecordsCount() {
    return mNumDeduped;
  }

  @Override
  public long getTotalNucleotides() {
    return mRecordWrapped.getTotalNucleotides();
  }

  @Override
  public long getFilteredRecordsCount() {
    return mRecordWrapped.getFilteredRecordsCount();
  }

  @Override
  public long getInvalidRecordsCount() {
    return mRecordWrapped.getInvalidRecordsCount();
  }

  @Override
  public long getOutputRecordsCount() {
    return getTotalRecordsCount() - getInvalidRecordsCount() - getDuplicateRecordsCount() - getFilteredRecordsCount();
  }

  @Override
  public long getTotalRecordsCount() {
    return mRecordWrapped.getTotalRecordsCount();
  }

  @Override
  public SAMFileHeader header() {
    return mRecordWrapped.header();
  }

  @Override
  public void close() throws IOException {
    mRecordWrapped.close();
  }

}
