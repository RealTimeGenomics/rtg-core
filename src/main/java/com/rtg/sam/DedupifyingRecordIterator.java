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
  public long getOverCoverageRecordsCount() {
    return mRecordWrapped.getOverCoverageRecordsCount();
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
    return mRecordWrapped.getOutputRecordsCount() - getDuplicateRecordsCount();
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
