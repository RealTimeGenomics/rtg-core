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
package com.rtg.variant;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.sam.ForcedInitialiser;
import com.rtg.util.Populator;

import htsjdk.samtools.SAMRecord;

/**
 * Populator that just passes through SAM records.
 */
public class SamRecordPopulator implements Populator<SAMRecord> {

  private static final SAMRecord OVERFLOW = new SAMRecord(null);

  private final ForcedInitialiser mRecordInitialiser;

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
    this(new DefaultForcedInitialiser());
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

