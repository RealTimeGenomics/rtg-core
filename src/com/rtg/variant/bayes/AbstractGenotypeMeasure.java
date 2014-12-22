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
package com.rtg.variant.bayes;

import com.reeltwo.jumble.annotations.TestClass;

/**
 */
@TestClass("com.rtg.variant.bayes.ArrayGenotypeMeasureTest")
public abstract class AbstractGenotypeMeasure implements GenotypeMeasure {
  private Integer mBest = null;
  private Double mBestPosterior = null;
  private Double mNonIdentity = null;
  private Hypotheses<?> mHyp;

  /**
   * Constructor
   * @param hyp the hypotheses this measure covers
   */
  public AbstractGenotypeMeasure(Hypotheses<?> hyp) {
    mHyp = hyp;
  }

  private void calculatePosteriors() {
    if (mBestPosterior != null) {
      return;
    }
    final int best = best();
    if (best == -1) {
      mBestPosterior = arithmetic().zero();
      mNonIdentity = arithmetic().zero();
      return;
    }
    double others = arithmetic().zero();
    final int reference = mHyp.reference();
    for (int i = 0; i < size(); i++) {
      if (i != best && i != reference) {
        others = arithmetic().add(others, measure(i));
      }
    }
    final double refMeasure = reference == -1 ? arithmetic().zero() : measure(reference);

    final double bestOther = best == reference ? others : arithmetic().add(others, refMeasure);
    mBestPosterior = arithmetic().divide(measure(best), bestOther);

    if (reference == -1) {
      mNonIdentity = arithmetic().zero();
    } else {
      final double refOther = best == reference ? others : arithmetic().add(others, measure(best));
      mNonIdentity = arithmetic().divide(refOther, refMeasure);
    }
  }


  @Override
  public double bestPosterior() {
    calculatePosteriors();
    return mBestPosterior;
  }

  @Override
  public double nonIdentityPosterior() {
    calculatePosteriors();
    return mNonIdentity;
  }

  private int findBest() {
    if (size() < 1) {
      return -1;
    }
    int best = 0;
    for (int i = 0; i < size(); i++) {
      if (arithmetic().gt(measure(i), measure(best))) {
        best = i;
      }
    }
    return best;
  }

  @Override
  public int best() {
    if (mBest == null) {
      mBest = findBest();
    }
    return mBest;
  }

  @Override
  public int reference() {
    return mHyp.reference();
  }

  @Override
  public Hypotheses<?> hypotheses() {
    return mHyp;
  }


  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append("AbstractGenotypeMeasure{").append("mBest=").append(mBest).append(" [");
    String join = "";
    for (int i = 0; i < size(); i++) {
      sb.append(join).append(measure(i));
      join = ", ";
    }
    sb.append("]}");
    return sb.toString();
  }
}
