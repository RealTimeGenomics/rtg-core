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

package com.rtg.variant.realign;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.util.VariantUtils;

/**
 * Probability function based on a q value and a reference nucleotide. It is
 * assumed that we are dealing with 4 nucleotides (0=A, ... 3=T) outside this
 * range returns a probability of 0.0 (Double.NEGATIVE
 *
 */
public class ProbabilityArray extends IntegralAbstract {

  private final double[] mProb;

  /**
   * @param probLn array of log probabilities.
   */
  public ProbabilityArray(final double[] probLn) {
    assert probLn.length == 4;
    //renormalize the probabilities
    double sum = Double.NEGATIVE_INFINITY;
    for (final double probLni : probLn) {
      sum = VariantUtils.logSumApproximation(sum, probLni);
    }
    mProb = new double[probLn.length];
    if (sum == Double.NEGATIVE_INFINITY) {
      for (int i = 0; i < probLn.length; ++i) {
        mProb[i] = 0.25;
      }
    } else {
      for (int i = 0; i < probLn.length; ++i) {
        final double d = probLn[i] - sum;
        mProb[i] = Math.exp(d);
      }
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertEquals(4, mProb.length);
    double sum = 0.0;
    for (double aMProb : mProb) {
      sum += aMProb;
    }
    Exam.assertEquals(1.0, sum, 0.0001);
    return true;
  }

  @Override
  public void toString(final StringBuilder sb) {
    //sb.append("ProbabilityArray");
    sb.append("[");
    for (int i = 0; i < mProb.length; ++i) {
      if (i > 0) {
        sb.append(" ");
      }
      sb.append(Utils.realFormat(mProb[i], 3));
    }
    sb.append("]");
  }
}
