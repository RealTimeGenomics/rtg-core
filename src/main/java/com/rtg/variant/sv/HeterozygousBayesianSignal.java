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

package com.rtg.variant.sv;


/**
 * Produces a bayesian signal for normal mappings (no structural variant).
 *
 */
public class HeterozygousBayesianSignal extends BayesianSignal {

  private final BayesianSignal mSignal1;
  private final BayesianSignal mSignal2;

  /**
   * Creates a new <code>HeterozygousBayesianSignal</code>. This
   * signal combines the distributions from two homozygous type
   * signals.
   *
   * @param signal1 a signal which will have it's distributions interpreted as a haploid distributions.
   * @param signal2 a signal which will have it's distributions interpreted as a haploid distributions.
   */
  public HeterozygousBayesianSignal(BayesianSignal signal1, BayesianSignal signal2) {
    mSignal1 = signal1;
    mSignal2 = signal2;
  }

  Distribution combineDistributions(Distribution a, Distribution b) {
    final int lo = a.lo();
    final int hi = a.hi();
    final int diameter = hi - lo;
    assert lo == b.lo();
    assert hi == b.hi() : hi + ":" + b.hi();
    final double[] distr = new double[diameter];
    for (int i = lo; i < hi; ++i) {
      distr[i - lo] = (a.get(i) + b.get(i)) / 2;
    }
    return new DistributionArray(lo, distr);
  }

  @Override
  Distribution leftArmProper(ReadGroupStats stats, boolean reverse) {
    return combineDistributions(mSignal1.leftArmProper(stats, reverse), mSignal2.leftArmProper(stats, reverse));
  }

  @Override
  Distribution leftArmDiscordant(ReadGroupStats stats, boolean reverse) {
    return combineDistributions(mSignal1.leftArmDiscordant(stats, reverse), mSignal2.leftArmDiscordant(stats, reverse));
  }
  @Override
  Distribution leftArmUnmated(ReadGroupStats stats, boolean reverse) {
    return combineDistributions(mSignal1.leftArmUnmated(stats, reverse), mSignal2.leftArmUnmated(stats, reverse));
  }
  @Override
  Distribution rightArmProper(ReadGroupStats stats, boolean reverse) {
    return combineDistributions(mSignal1.rightArmProper(stats, reverse), mSignal2.rightArmProper(stats, reverse));
  }
  @Override
  Distribution rightArmDiscordant(ReadGroupStats stats, boolean reverse) {
    return combineDistributions(mSignal1.rightArmDiscordant(stats, reverse), mSignal2.rightArmDiscordant(stats, reverse));
  }
  @Override
  Distribution rightArmUnmated(ReadGroupStats stats, boolean reverse) {
    return combineDistributions(mSignal1.rightArmUnmated(stats, reverse), mSignal2.rightArmUnmated(stats, reverse));
  }
}
