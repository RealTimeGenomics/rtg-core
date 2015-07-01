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
 * Maintains a set of counts over the possible nucleotides.
 */
public class AlleleStatisticsInt extends AlleleStatistics<AlleleStatisticsInt> {

  /**
   * Counts of different nucleotides.
   */
  private final int[] mCountsForwards;
  private final int[] mCountsBackwards;
  private final int[] mCountsMated;
  private final int[] mCountsUnmated;
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
  public AlleleStatisticsInt(final Description description) {
    super(description);
    mCountsForwards = new int[description.size()];
    mCountsBackwards = new int[description.size()];
    mCountsMated = new int[description.size()];
    mCountsUnmated = new int[description.size()];
    //mCountsSingle = new int[description.size()];
    mErrors = new double[description.size()];
    //Arrays.fill(mCountsForwards, 0);
    //Arrays.fill(mCountsBackwards, 0);
    //Arrays.fill(mErrors, 0.0);
  }

  /**
   * Increment the count for index making allowance for errors in read extraction.
   * @param distribution values to increment
   * @param bestHyp the code for the allele that best represents the read.
   * @param e error for current read
   */
  public void increment(final EvidenceInterface distribution, int bestHyp, double e) {
    assert getDescription().valid(bestHyp) : bestHyp;
    //System.err.println("r=" + r + " q=" + q + " e=" + e);
    //System.err.println("Counts increment index=" + index + " e=" + e);

    if (distribution.isForward()) {
      mCountsForwards[bestHyp]++;
    } else {
      mCountsBackwards[bestHyp]++;
    }
    if (distribution.isReadPaired()) {
      if (distribution.isMated()) {
        mCountsMated[bestHyp]++;
      } else {
        mCountsUnmated[bestHyp]++;
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
  public double count(final int index) {
    return mCountsForwards[index] + mCountsBackwards[index];
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
  Double strandBias(int allele) {
    final long trials = MathUtils.round(count(allele));
    final long observed = mCountsForwards[allele];
    return MathUtils.hoeffdingPhred(trials, observed, 0.5);
  }

  @Override
  Double unmatedBias(int allele, double unmatedProbability) {
    final int trials = mCountsMated[allele] + mCountsUnmated[allele];
    final int observed = mCountsUnmated[allele];
    return MathUtils.hoeffdingPhred(trials, observed, unmatedProbability);
  }


  @Override
  public AlleleStatisticsInt remap(Description newDescription, int[] mapping) {
    final AlleleStatisticsInt newCounts = new AlleleStatisticsInt(newDescription);
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

//  String debug() {
//    final StringBuilder sb = new StringBuilder();
//    sb.append("f|b|n!");
//    for (int i = 0; i < mDescription.size(); i++) {
//      sb.append("").append(mCountsForwards[i]).append("|").append(mCountsBackwards[i]).append("|").append(mDescription.name(i)).append("!");
//    }
//    final int lastIndex = mCountsForwards.length - 1;
//    sb.append("").append(mCountsForwards[lastIndex]).append("|").append(mCountsBackwards[lastIndex]).append("|").append("other").append("!");
//    return sb.toString();
//  }
}
