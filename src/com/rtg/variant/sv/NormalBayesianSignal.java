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
