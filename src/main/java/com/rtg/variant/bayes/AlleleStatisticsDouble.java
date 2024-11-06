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

import com.rtg.util.MathUtils;

/**
 * Provides statistics per allele at double resolution for coverage.
 */
public class AlleleStatisticsDouble extends AlleleStatistics<AlleleStatisticsDouble> {

  /**
   * Counts of different nucleotides.
   */
  private final double[] mCountsForwards1;
  private final double[] mCountsForwards2;
  private final double[] mCountsBackwards1;
  private final double[] mCountsBackwards2;
  private final double[] mCountsMated;
  private final double[] mCountsUnmated;

  /**
   * Sum of error probability for each nucleotide.
   */
  private final double[] mErrors;

  /**
   * Product of base quality errors.
   */
  private final double[] mQualityProduct;


  /**
   * Create a new empty count object.
   * @param description set of possible values
   */
  public AlleleStatisticsDouble(final Description description) {
    super(description);
    mCountsForwards1 = new double[description.size()];
    mCountsForwards2 = new double[description.size()];
    mCountsBackwards1 = new double[description.size()];
    mCountsBackwards2 = new double[description.size()];
    mCountsMated = new double[description.size()];
    mCountsUnmated = new double[description.size()];
    mErrors = new double[description.size()];
    mQualityProduct = new double[description.size()];
  }


  /**
   * Increment the count for index making allowance for errors in read extraction.
   * @param distribution values to increment
   * @param bestHyp the code for the allele that best represents the read.
   * @param e error for current read
   * @param coverage amount of coverage for this increment
   */
  public void increment(final EvidenceInterface distribution, int bestHyp, double e, double coverage) {
    assert getDescription().valid(bestHyp) : bestHyp;
    if (distribution.isForward()) {
      if (distribution.isFirst()) {
        mCountsForwards1[bestHyp] += coverage;
      } else {
        mCountsForwards2[bestHyp] += coverage;
      }
    } else {
      if (distribution.isFirst()) {
        mCountsBackwards1[bestHyp] += coverage;
      } else {
        mCountsBackwards2[bestHyp] += coverage;
      }
    }
    if (distribution.isReadPaired()) {
      if (distribution.isMated()) {
        mCountsMated[bestHyp] += coverage;
      } else {
        mCountsUnmated[bestHyp] += coverage;
      }
    }
    mErrors[bestHyp] += e;
    mQualityProduct[bestHyp] += MathUtils.phred(e);
  }

  @Override
  public double count(final int index) {
    return forward(index) + backward(index);
  }

  @Override
  public double forward1(final int index) {
    return mCountsForwards1[index];
  }

  @Override
  public double forward2(final int index) {
    return mCountsForwards2[index];
  }

  @Override
  public double backward1(final int index) {
    return mCountsBackwards1[index];
  }

  @Override
  public double backward2(final int index) {
    return mCountsBackwards2[index];
  }

  @Override
  public double error(final int index) {
    return mErrors[index];
  }

  @Override
  public double qa(final int index) {
    return mQualityProduct[index];
  }

  @Override
  Double strandBias(int allele) {
    final double trials = count(allele);
    final double observed = forward(allele);
    return MathUtils.hoeffdingPhred(trials, observed, 0.5);
  }

  @Override
  Double unmatedBias(int allele, double unmatedProbability) {
    final double trials = mCountsMated[allele] + mCountsUnmated[allele];
    final double observed = mCountsUnmated[allele];
    return MathUtils.hoeffdingPhred(trials, observed, unmatedProbability);
  }

  @Override
  public AlleleStatisticsDouble remap(Description newDescription, int[] mapping) {
    final AlleleStatisticsDouble newCounts = new AlleleStatisticsDouble(newDescription);
    for (int oldI = 0; oldI < mapping.length; ++oldI) {
      final int newI = mapping[oldI];
      if (newI >= 0) {
        newCounts.mCountsForwards1[newI] += mCountsForwards1[oldI];
        newCounts.mCountsForwards2[newI] += mCountsForwards2[oldI];
        newCounts.mCountsBackwards1[newI] += mCountsBackwards1[oldI];
        newCounts.mCountsBackwards2[newI] += mCountsBackwards2[oldI];
        newCounts.mCountsMated[newI] += mCountsMated[oldI];
        newCounts.mCountsUnmated[newI] += mCountsUnmated[oldI];
        newCounts.mErrors[newI] += mErrors[oldI];
        newCounts.mQualityProduct[newI] += mQualityProduct[oldI];
      }
    }
    return newCounts;
  }
}
