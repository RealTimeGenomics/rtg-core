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
}
