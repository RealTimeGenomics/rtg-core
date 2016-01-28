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
 * Approximate probability of allele balance with Hoeffding's inequality
 * @author kurt
 */
public class HoeffdingAlleleBalance extends AbstractAlleleBalance {
  /**
   * @param expected expected allele balance
   */
  public HoeffdingAlleleBalance(double expected) {
    super(expected);
  }

  private double alleleBalanceLn(double p, double trials, double count) {
    return MathUtils.hoeffdingLn(trials, count, p);
  }

  @Override
  double alleleBalanceHeterozygousLn(double p, double trials, double observed, double observedAlt) {
    return alleleBalanceLn(p, trials, observed) + alleleBalanceLn(p, trials, observedAlt);
  }
}
