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

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Produces a bayesian signal for normal mappings (no structural variant).
 *
 */
@TestClass(value = {"com.rtg.variant.sv.SvToolTaskTest"})
public class NormalBayesianSignal extends BayesianSignal {

  private final double mCoverage;

  /**
   * Creates a new <code>NormalBayesianSignal</code>.
   *
   * @param coverage the level of coverage that this signal detects as being normal.
   */
  public NormalBayesianSignal(int coverage) {
    mCoverage = coverage / 2.0;
  }

  @Override
  Distribution leftArmProper(ReadGroupStats stats, boolean reverse) {
    final int lo = stats.lo();
    final int hi = stats.hi();
    return new DistributionConstant(lo, hi, mCoverage * stats.properRate() + stats.properRandomRate());
  }
  @Override
  Distribution leftArmDiscordant(ReadGroupStats stats, boolean reverse) {
    final int lo = stats.lo();
    final int hi = stats.hi();
    return new DistributionConstant(lo, hi, stats.discordantRate());
  }
  @Override
  Distribution leftArmUnmated(ReadGroupStats stats, boolean reverse) {
    final int lo = stats.lo();
    final int hi = stats.hi();
    return new DistributionConstant(lo, hi, stats.unmatedRate());
  }
  @Override
  Distribution rightArmProper(ReadGroupStats stats, boolean reverse) {
    return leftArmProper(stats, reverse);
  }
  @Override
  Distribution rightArmDiscordant(ReadGroupStats stats, boolean reverse) {
    return leftArmDiscordant(stats, reverse);
  }
  @Override
  Distribution rightArmUnmated(ReadGroupStats stats, boolean reverse) {
    return leftArmUnmated(stats, reverse);
  }
}
