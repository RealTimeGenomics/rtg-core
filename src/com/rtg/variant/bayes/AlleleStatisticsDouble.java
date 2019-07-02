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
    final long trials = MathUtils.round(count(allele));
    final long observed = MathUtils.round(forward(allele));
    return MathUtils.hoeffdingPhred(trials, observed, 0.5);
  }

  @Override
  Double unmatedBias(int allele, double unmatedProbability) {
    final long trials = MathUtils.round(mCountsMated[allele] + mCountsUnmated[allele]);
    final long observed = MathUtils.round(mCountsUnmated[allele]);
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
