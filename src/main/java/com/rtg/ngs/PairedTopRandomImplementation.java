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
