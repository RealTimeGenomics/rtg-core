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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.reference.Ploidy;
import com.rtg.variant.util.arithmetic.LogPossibility;

/**
 */
@TestClass({"com.rtg.variant.bayes.BinomialAlleleBalanceTest", "com.rtg.variant.bayes.HoeffdingAlleleBalanceTest"})
public abstract class AbstractAlleleBalance implements AlleleBalanceProbability {
  final double mExpected;

  /**
   * @param expected the expected allele balance
   */
  protected AbstractAlleleBalance(double expected) {
    this.mExpected = expected;
  }

  @Override
  public double alleleBalanceLn(int i, Hypotheses<?> hypotheses, Statistics<?> statistics) {
    if (hypotheses.ploidy() != Ploidy.DIPLOID
      && hypotheses.ploidy() != Ploidy.HAPLOID
      ) {
      return LogPossibility.SINGLETON.one();
    }
    final double trials = statistics instanceof StatisticsDouble
      ? ((StatisticsDouble) statistics).exactCoverage()
      : statistics.coverage();

    if (trials == 0) {
      return LogPossibility.SINGLETON.one();
    }
    final int a = hypotheses.code().a(i);
    final int b = hypotheses.code().bc(i);
    final AlleleStatistics<?> counts = statistics.counts();
    final double vac = counts.count(a) - counts.error(a);
    if (a == b) {
       //double error = statistics.totalError() / statistics.coverage();
       //alleleBalanceHomozygousLn(1.0 - error, trials, vac);
      return LogPossibility.SINGLETON.one();
    }
    final double bCount = counts.count(b) - counts.error(b);
    return alleleBalanceHeterozygousLn(mExpected, trials, vac, bCount);
  }

  abstract double alleleBalanceHeterozygousLn(double p, double trials, double observed, double observedAlt);

}
