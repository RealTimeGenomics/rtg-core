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
