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

import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Mock Genotype measure for testing
 */
public class MockGenotypeMeasure implements GenotypeMeasure {
  final int mReference;
  final int mBest;
  final double mPosterior;
  final double mNonIdentityPosterior;

  /**
   * Mock measure containing just a bestPosterior
   * @param posterior the posterior to use
   */
  public MockGenotypeMeasure(double posterior) {
    this(0, 0, posterior, 0);
  }

  /**
   *
   * @param reference reference call index
   * @param best best call index
   * @param posterior best call posterior
   * @param nonIdentityPosterior posterior not equal to reference
   */
  public MockGenotypeMeasure(int reference, int best, double posterior, double nonIdentityPosterior) {
    mReference = reference;
    mBest = best;
    mPosterior = posterior;
    mNonIdentityPosterior = nonIdentityPosterior;

  }
  @Override
  public PossibilityArithmetic arithmetic() {
    return LogApproximatePossibility.SINGLETON;
  }

  @Override
  public double measure(int hypothesis) {
    return LogApproximatePossibility.SINGLETON.prob2Poss(1.0 / 10.0);
  }

  @Override
  public int size() {
    throw new UnsupportedOperationException();
  }

  @Override
  public double bestPosterior() {
    return mPosterior;
  }

  @Override
  public double nonIdentityPosterior() {
    return mNonIdentityPosterior;
  }

  @Override
  public int best() {
    return mBest;
  }

  @Override
  public int reference() {
    return mReference;
  }

  @Override
  public Hypotheses<?> hypotheses() {
    return new HypothesesSnp(arithmetic(), GenomePriorParams.builder().create(), false, mReference);
  }
}
