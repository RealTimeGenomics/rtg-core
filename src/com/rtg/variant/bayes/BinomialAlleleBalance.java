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
 * @author kurt
 */
public class BinomialAlleleBalance extends AbstractAlleleBalance {

  private static double logBinomial(final double p, final double nn, final double n) {
    assert p >= 0.0 && p <= 1.0;
    assert n >= 0;
    assert nn >= 0;
    final double m = nn - n;
    final double res = n * Math.log(p) + m * Math.log(1.0 - p) + ChiSquared.lgamma(nn + 1) - ChiSquared.lgamma(n + 1) - ChiSquared.lgamma(m + 1);
    assert res <= 0 && !Double.isNaN(res);
    return res;
  }

  @Override
  double alleleBalanceHeterozygousLn(double p, double trials, double count, double countAlt) {
    return alleleBalanceLn(p, trials, count) + alleleBalanceLn(p, trials, countAlt);
  }

  private double alleleBalanceLn(double p, double trials, double count) {
    return logBinomial(p, trials, count);
  }
}
