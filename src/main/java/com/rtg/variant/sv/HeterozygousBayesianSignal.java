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
