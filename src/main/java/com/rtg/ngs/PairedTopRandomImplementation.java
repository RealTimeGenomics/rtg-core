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
package com.rtg.ngs;

import java.util.Random;

/**
 *
 */
public class PairedTopRandomImplementation {

  private HitRecord[] mRecords;

  private final Random mRandom;

  private final int mNumReads;

  /*
   * We do not want to add any unnecessary fields into this class, as there will
   * be NUM_READS of them, and every byte will cause considerable memory usage
   */
  static final class HitRecord {
    int mTemplateId;
    //int mReadId;
    int mZeroBasedTemplateStart;
    boolean mReverse;

    int mZeroBasedMateTemplateStart;
    boolean mMateReverse;
    int mComboScore;
//    int mInferredInsertSize;
    //maybe add score here to make re-alignment more efficient

    private HitRecord() { }

  }

  /**
   * @param numReads number of read pairs
   */
  public PairedTopRandomImplementation(int numReads) {
    mNumReads = numReads;
    mRandom = new Random();
  }

  PairedTopRandomImplementation(int numReads, long seed) {
    mNumReads = numReads;
    mRandom = new Random(seed);
  }

  /**
   * Initialise the Top Random implementation
   */
  public void initialize() {
    mRecords = new HitRecord[mNumReads];
  }

  /**
   * Update top random structure with a <code>1/n</code> probability.
   * @param readId the read id
   * @param templateId the template id of the sequence this read maps to.
   * @param zeroBasedTemplateStart start position <code>(0-based)</code> of first (sequencing wise) arm
   * @param reverse true if first arm was mapped reverse complement
   * @param zeroBasedMateTemplateStart start position <code>(0-based)</code> of second (sequencing wise) arm
   * @param mateReverse true if second arm was mapped reverse complement
   * @param comboScore the score of this combo. Used to avoid race conditions.
   * @param n number of records seen at alignments current combo score
   */
  public void update(int readId, int templateId, int zeroBasedTemplateStart, boolean reverse, int zeroBasedMateTemplateStart, boolean mateReverse, int comboScore, int n) {
    //protect against race condition
    if (mRecords[readId] != null && comboScore > mRecords[readId].mComboScore) {
      return;
    }
    if (n != 1) {
      if (mRandom.nextInt(n) != 0) {
        return;
      }
    }
    final HitRecord rec = getRecord(readId);
    //rec.mReadId = readId;
    rec.mReverse = reverse;
    rec.mTemplateId = templateId;
    rec.mZeroBasedTemplateStart = zeroBasedTemplateStart;
    rec.mZeroBasedMateTemplateStart = zeroBasedMateTemplateStart;
    rec.mMateReverse = mateReverse;
    rec.mComboScore = comboScore;
  }

  private HitRecord getRecord(int index) {
    if (mRecords[index] == null) {
      mRecords[index] = new HitRecord();
    }
    return mRecords[index];
  }

  HitRecord[] getRecords() {
    return mRecords;
  }

  void finish() {
    mRecords = null;
  }
}
