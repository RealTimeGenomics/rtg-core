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
  private final int[] mCountsForwards1;
  private final int[] mCountsForwards2;
  private final int[] mCountsBackwards1;
  private final int[] mCountsBackwards2;
  private final int[] mCountsMated;
  private final int[] mCountsUnmated;

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
  public AlleleStatisticsInt(final Description description) {
    super(description);
    mCountsForwards1 = new int[description.size()];
    mCountsForwards2 = new int[description.size()];
    mCountsBackwards1 = new int[description.size()];
    mCountsBackwards2 = new int[description.size()];
    mCountsMated = new int[description.size()];
    mCountsUnmated = new int[description.size()];
    mErrors = new double[description.size()];
    mQualityProduct = new double[description.size()];
  }

  /**
   * Increment the count for index making allowance for errors in read extraction.
   * @param distribution values to increment
   * @param bestHyp the code for the allele that best represents the read.
   * @param e error for current read
   */
  public void increment(final EvidenceInterface distribution, int bestHyp, double e) {
    assert getDescription().valid(bestHyp) : bestHyp;

    if (distribution.isForward()) {
      if (distribution.isFirst()) {
        ++mCountsForwards1[bestHyp];
      } else {
        ++mCountsForwards2[bestHyp];
      }
    } else {
      if (distribution.isFirst()) {
        ++mCountsBackwards1[bestHyp];
      } else {
        ++mCountsBackwards2[bestHyp];
      }
    }
    if (distribution.isReadPaired()) {
      if (distribution.isMated()) {
        ++mCountsMated[bestHyp];
      } else {
        ++mCountsUnmated[bestHyp];
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
    final double observed = forward(allele);
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

//  String debug() {
//    final StringBuilder sb = new StringBuilder();
//    sb.append("f|b|n!");
//    for (int i = 0; i < mDescription.size(); ++i) {
//      sb.append("").append(mCountsForwards[i]).append("|").append(mCountsBackwards[i]).append("|").append(mDescription.name(i)).append("!");
//    }
//    final int lastIndex = mCountsForwards.length - 1;
//    sb.append("").append(mCountsForwards[lastIndex]).append("|").append(mCountsBackwards[lastIndex]).append("|").append("other").append("!");
//    return sb.toString();
//  }
}
