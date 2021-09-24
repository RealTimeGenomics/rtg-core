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
package com.rtg.variant.coverage;

import com.rtg.sam.MateInfo;
import com.rtg.sam.ReaderRecord;
import com.rtg.util.intervals.SequenceIdLocusSimple;

import htsjdk.samtools.SAMRecord;

/**
 * Handles the standard mate start/sequence information as well as the record chaining
 */
public abstract class AbstractMateInfoReaderRecord<T extends AbstractMateInfoReaderRecord<T>> extends SequenceIdLocusSimple implements MateInfo, ReaderRecord<T> {
  protected final boolean mMated;
  protected final int mMateSequenceId;
  protected final int mFragmentLength;
  protected final int mGenome;
  private T mNextInChain;

  /**
   * Constructor
   * @param sam the sam record this is based on
   * @param genome the genome id
   */
  public AbstractMateInfoReaderRecord(SAMRecord sam, int genome) {
    super(sam.getReferenceIndex(), sam.getAlignmentStart() - 1, sam.getAlignmentEnd());
    mMated = sam.getReadPairedFlag() && sam.getProperPairFlag();
    if (mMated) {
      mFragmentLength = sam.getInferredInsertSize();
      mMateSequenceId = sam.getMateReferenceIndex();
    } else {
      mMateSequenceId = -1;
      mFragmentLength = -1;
    }
    mGenome = genome;
  }

  @Override
  public T chain() {
    return mNextInChain;
  }

  @Override
  public int getGenome() {
    return mGenome;
  }

  @Override
  public void setNextInChain(T rec) {
    mNextInChain = rec;
  }

  @Override
  public int getMateSequenceId() {
    return mMateSequenceId;
  }

  @Override
  public int getFragmentLength() {
    return mFragmentLength;
  }

  @Override
  public boolean isMated() {
    return mMated;
  }
}
