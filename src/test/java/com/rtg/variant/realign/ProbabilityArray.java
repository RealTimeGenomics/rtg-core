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
