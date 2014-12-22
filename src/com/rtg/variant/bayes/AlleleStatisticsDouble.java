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

import java.util.Arrays;

import com.rtg.util.MathUtils;

/**
 * Provides statistics per allele at double resolution for coverage
 */
public class AlleleStatisticsDouble extends AlleleStatistics<AlleleStatisticsDouble> {

  /**
   * Counts of different nucleotides.
   */
  private final double[] mCountsForwards;
  private final double[] mCountsBackwards;
  private final double[] mCountsMated;
  private final double[] mCountsUnmated;
  //private final int[] mCountsSingle;
  //private static final double LOG10_E = Math.log(Math.E);


  /**
   * Sum of error probability for each nucleotide.
   */
  private final double[] mErrors;


  /**
   * Create a new empty count object.
   * @param description set of possible values
   */
  public AlleleStatisticsDouble(final Description description) {
    super(description);
    mCountsForwards = new double[description.size()];
    mCountsBackwards = new double[description.size()];
    mCountsMated = new double[description.size()];
    mCountsUnmated = new double[description.size()];
    //mCountsSingle = new int[description.size()];
    mErrors = new double[description.size()];
    Arrays.fill(mCountsForwards, 0);
    Arrays.fill(mCountsBackwards, 0);
    Arrays.fill(mErrors, 0.0);
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
    //System.err.println("r=" + r + " q=" + q + " e=" + e);
    //System.err.println("Counts increment index=" + index + " e=" + e);

    if (distribution.isForward()) {
      mCountsForwards[bestHyp] += coverage;
    } else {
      mCountsBackwards[bestHyp] += coverage;
    }
    if (distribution.isReadPaired()) {
      if (distribution.isMated()) {
        mCountsMated[bestHyp] += coverage;
      } else {
        mCountsUnmated[bestHyp] += coverage;
      }
    }
    mErrors[bestHyp] += e;
  }

  /**
   * Get the current count for the specified index.
   * @param index whose value to get.
   * @return the count.
   */
  @Override
  public int count(final int index) {
    return (int) MathUtils.round(mCountsForwards[index] + mCountsBackwards[index]);
  }

  /**
   * Get the current accumulated for the specified index.
   * @param index whose value to get.
   * @return the accumulated error.
   */
  @Override
  public double error(final int index) {
    return mErrors[index];
  }

  @Override
  Double alleleBalanceHomozygous(int allele, int total) {
    final int observed = (int) MathUtils.round(mCountsForwards[allele] + mCountsBackwards[allele]);
    return MathUtils.hoeffdingPhred(total, observed, 1.0);
  }

  @Override
  Double alleleBalance(int allele1, int allele2) {
    final int trials = (int) MathUtils.round(mCountsForwards[allele1] + mCountsBackwards[allele2] + mCountsForwards[allele2] + mCountsBackwards[allele1]);
    final int observed = (int) MathUtils.round(mCountsForwards[allele1] + mCountsBackwards[allele1]);
    return MathUtils.hoeffdingPhred(trials, observed, 0.5);
  }

  @Override
  Double strandBias(int allele) {
    final int trials = (int) MathUtils.round(mCountsForwards[allele] + mCountsBackwards[allele]);
    final int observed = (int) MathUtils.round(mCountsForwards[allele]);
    return MathUtils.hoeffdingPhred(trials, observed, 0.5);
  }

  @Override
  Double unmatedBias(int allele, double unmatedProbability) {
    final int trials = (int) MathUtils.round(mCountsMated[allele] + mCountsUnmated[allele]);
    final int observed = (int) MathUtils.round(mCountsUnmated[allele]);
    return MathUtils.hoeffdingPhred(trials, observed, unmatedProbability);
  }

  @Override
  public AlleleStatisticsDouble remap(Description newDescription, int[] mapping) {
    final AlleleStatisticsDouble newCounts = new AlleleStatisticsDouble(newDescription);
    for (int oldI = 0; oldI < mapping.length; oldI++) {
      final int newI = mapping[oldI];
      newCounts.mCountsForwards[newI] += mCountsForwards[oldI];
      newCounts.mCountsBackwards[newI] += mCountsBackwards[oldI];
      newCounts.mCountsMated[newI] += mCountsMated[oldI];
      newCounts.mCountsUnmated[newI] += mCountsUnmated[oldI];
      newCounts.mErrors[newI] += mErrors[oldI];
    }
    return newCounts;
  }

}
