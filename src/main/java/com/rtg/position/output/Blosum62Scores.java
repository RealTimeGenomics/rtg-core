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
package com.rtg.position.output;

import com.rtg.mode.Frame;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Pseudo probabilities when using Blosum62 scoring.
 * Returns 1.0 for gaps that are ok 0.0 otherwise.
 * The consequent values are computed from this.
 */
public class Blosum62Scores extends IntegralAbstract implements GapScorer {

  private final int mMaxGap;

  /** Gaps less or equal to this are permitted to have indels. */
  private final int mIGap;

  /**
   * @param maxGap maximum allowed gap.
   */
  public Blosum62Scores(final int maxGap) {
    mMaxGap = maxGap;
    mIGap = mMaxGap < 10 ? 0 : mMaxGap - 10; //TODO make this real
  }

  private boolean allowed(final int i, final int j) {
    return i == j && i <= mMaxGap || i <= mIGap && j <= mIGap;
  }


  @Override
  public double scoreMax(int buildStart, int buildEnd, int queryStart, int queryEnd) {
    final int i = buildEnd - buildStart;
    return i <= mMaxGap ? scoreInner(i, i) : Double.NEGATIVE_INFINITY;
  }

  private double scoreInner(final int i, final int j) {
    assert i >= 0 && j >= 0;
    if (!allowed(i, j)) {
      return Double.NEGATIVE_INFINITY;
    }
    if (i == j) {
      return i * Blosum62.EXPECTED;
    }
    final int min;
    final int del;
    if (i < j) {
      min = i;
      del = j - i;
    } else {
      min = j;
      del = i - j;
    }
    assert del >= 1;
    final double sc = min * Blosum62.EXPECTED + Blosum62.GAP + del * Blosum62.EXTEND;
    assert sc < 0.0;
    return sc;
  }

  @Override
  public double score(int buildStart, int buildEnd, int queryStart, int queryEnd) {
    final int i = buildEnd - buildStart;
    final int j = queryEnd - queryStart;
    //System.err.println("i=" + i + " j=" + j);
    if (i < 0) {
      //overlap
      return i == j ? 0.0 : Double.NEGATIVE_INFINITY;
    }
    if (j < 0) {
      return Double.NEGATIVE_INFINITY;
    }
    return scoreInner(i, j);
  }

  @Override
  public int maxDelta() {
    return mIGap;
  }

  @Override
  public int minDelta() {
    return -mIGap;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Blosum62Probabilities");
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mMaxGap >= 1);
    Exam.assertTrue(mIGap >= 0);
    return true;
  }

  @Override
  public void setInternalBuildId(int internalId) {
  }

  @Override
  public void setQueryId(int internalId, Frame frame) {
  }


  @Override
  public void close() {
  }

}
