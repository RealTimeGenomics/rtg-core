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
package com.rtg.variant.bayes;

import com.reeltwo.jumble.annotations.TestClass;

/**
 */
@TestClass("com.rtg.variant.bayes.ArrayGenotypeMeasureTest")
public abstract class AbstractGenotypeMeasure implements GenotypeMeasure {
  private Integer mBest = null;
  private Double mBestPosterior = null;
  private Double mNonIdentity = null;
  private final Hypotheses<?> mHyp;

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
    for (int i = 0; i < size(); ++i) {
      if (i != best && i != reference) {
        others = arithmetic().add(others, measure(i));
      }
    }
    final double refMeasure = reference == Hypotheses.NO_HYPOTHESIS ? arithmetic().zero() : measure(reference);

    final double bestOther = best == reference ? others : arithmetic().add(others, refMeasure);
    mBestPosterior = arithmetic().divide(measure(best), bestOther);

    if (reference == Hypotheses.NO_HYPOTHESIS) {
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
      return Hypotheses.NO_HYPOTHESIS;
    }
    int best = 0;
    for (int i = 0; i < size(); ++i) {
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
    for (int i = 0; i < size(); ++i) {
      sb.append(join).append(hypotheses().name(i)).append("=").append(measure(i));
      join = ", ";
    }
    sb.append("]}");
    return sb.toString();
  }
}
