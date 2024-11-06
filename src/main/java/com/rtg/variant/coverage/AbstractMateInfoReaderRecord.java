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
