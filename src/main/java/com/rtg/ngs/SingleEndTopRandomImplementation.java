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
public class SingleEndTopRandomImplementation {

  private HitRecord[] mRecords;

  private final Random mRandom;

  private final int mNumRecords;

  /*
   * We do not want to add any unnecessary fields into this class, as there will
   * be 2 * NUM_UNMATED_READS of them (paired-end), or NUM_READS (single-end),
   * and every byte will cause considerable memory usage
   */
  static final class HitRecord {
    int mTemplateId;
    //int mReadId;
    int mZeroBasedTemplateStart;
    boolean mReverse;
    int mAlignScore;

    private HitRecord() { }

  }

  /**
   * @param numRecords maximum number of records to store (
   */
  public SingleEndTopRandomImplementation(int numRecords) {
    mNumRecords = numRecords;
    mRandom = new Random();
  }

  SingleEndTopRandomImplementation(int numRecords, long seed) {
    mNumRecords = numRecords;
    mRandom = new Random(seed);
  }

  /**
   * Initialize the Top Random implementation
   */
  public void initialize() {
    mRecords = new HitRecord[mNumRecords];
  }

  /**
   * Update top random structure with a <code>1/n</code> probability.
   * @param readId the read id
   * @param templateId the template id of the sequence this read maps to.
   * @param zeroBasedTemplateStart start position <code>(0-based)</code> of first (sequencing wise) arm
   * @param reverse true if first arm was mapped reverse complement
   * @param alignScore alignment score of hit
   * @param n number of records seen at alignments current combo score
   */
  public void update(int readId, int templateId, int zeroBasedTemplateStart, boolean reverse, int alignScore, int n) {
    //protect against race condition
    if (mRecords[readId] != null && alignScore > mRecords[readId].mAlignScore) {
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
    rec.mAlignScore = alignScore;
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
