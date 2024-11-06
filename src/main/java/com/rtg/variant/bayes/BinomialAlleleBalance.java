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

import com.rtg.util.ChiSquared;

/**
 * Approximate probability of allele balance as a binomial trial.
 */
public class BinomialAlleleBalance extends AbstractAlleleBalance {

  /**
   * @param expected the expected allele balance
   */
  public BinomialAlleleBalance(double expected) {
    super(expected);
  }

  // C(total, count) * p^count * (1-p)^(total-count)
//  private static double logBinomial(final double p, final double total, final double count) {
//    assert p >= 0.0 && p <= 1.0;
//    assert count >= 0;
//    assert total >= 0;
//    final double m = total - count;
//    final double res = count * Math.log(p) + m * Math.log(1.0 - p) + ChiSquared.lgamma(total + 1) - ChiSquared.lgamma(count + 1) - ChiSquared.lgamma(m + 1);
//    assert res <= 0 && !Double.isNaN(res);
//    return res;
//  }

  // Compute a probability in log space using a binomial style distribution for the probability of seeing
  // count1/total having expected probability p and count2/total having expected probability 1-p
  private static double logBinomial(final double p, final double total, final double count1, final double count2) {
    assert p >= 0.0 && p <= 1.0;
    assert count1 >= 0;
    assert count2 >= 0;
    assert total >= 0;
    final double res = (count1 + total - count2) * Math.log(p) + (total - count1 + count2) * Math.log(1.0 - p)
      + 2 * ChiSquared.lgamma(total + 1)
      - ChiSquared.lgamma(count1 + 1) - ChiSquared.lgamma(total - count1 + 1)
      - ChiSquared.lgamma(count2 + 1) - ChiSquared.lgamma(total - count2 + 1);
    assert res <= 0 && !Double.isNaN(res);
    return res;
  }

  @Override
  double alleleBalanceHeterozygousLn(double p, double trials, double count, double countAlt) {
    // Compute bin(p, trials, count) * bin(1 - p, trials, countAlt)
    return logBinomial(p, trials, count, countAlt);
  }

  @Override
  public String toString() {
    return "binomial(" + mExpected + ")";
  }

}
