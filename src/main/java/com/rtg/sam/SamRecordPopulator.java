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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.Populator;

import htsjdk.samtools.SAMRecord;

/**
 * Populator that just passes through SAM records.
 */
public class SamRecordPopulator implements Populator<SAMRecord> {

  @TestClass(value = {"com.rtg.sam.ThreadedMultifileIteratorTest"})
  private static class DefaultForcedInitialiser implements ForcedInitialiser {
    @Override
    public void forceLazyInit(SAMRecord record) {
      record.getCigar();
      record.getCigarString();
      record.getAlignmentBlocks();
      record.getAlignmentEnd();
    }
  }

  /** Singleton for default record initialization */
  public static final ForcedInitialiser DEFAULT_INITIALIZER = new DefaultForcedInitialiser();

  private static final SAMRecord OVERFLOW = new SAMRecord(null);

  private final ForcedInitialiser mRecordInitialiser;

  /**
   * Populator for SAM records using an initializer. Can be given null to
   * skip initialization.
   * @param initializer the initializer
   */
  public SamRecordPopulator(final ForcedInitialiser initializer) {
    mRecordInitialiser = initializer;
  }

  /** Populator for SAM records with default initializer. */
  public SamRecordPopulator() {
    this(DEFAULT_INITIALIZER);
  }


  @Override
  public SAMRecord populate(final SAMRecord rec) {
    if (mRecordInitialiser != null) {
      mRecordInitialiser.forceLazyInit(rec);
    }
    return rec;
  }

  @Override
  public SAMRecord overflow(int position, int length) {
    return OVERFLOW;
  }
}

