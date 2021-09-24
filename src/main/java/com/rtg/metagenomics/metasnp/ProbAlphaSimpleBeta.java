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
package com.rtg.metagenomics.metasnp;

import java.util.Arrays;

/**
 */
public class ProbAlphaSimpleBeta implements ProbAlpha {
  private final double[] mBeta;

  /**
   * Alpha probability backed by a per strain variant prob
   * @param beta strain mutant probabilities in prob space
   */
  public ProbAlphaSimpleBeta(double[] beta) {
    mBeta = beta;
  }
  @Override
  public double pAlpha(int referenceAllele, int[] strainVariants) {
    double prob = 1;
    for (int i = 0; i < strainVariants.length; ++i) {
      final Integer assignment = strainVariants[i];
      prob *= assignment == referenceAllele ?  1 - mBeta[i] : mBeta[i];
    }
    return prob;
  }

  @Override
  public String toString() {
    return "ProbAlphaSimpleBeta{" + "beta: " + Arrays.toString(mBeta) + '}';
  }
}
