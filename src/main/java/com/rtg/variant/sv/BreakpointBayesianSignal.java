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
 * Produces a bayesian signal for a simple breakpoint.
 * This covers a duplication recipient.
 *
 */
public class BreakpointBayesianSignal extends BayesianSignal {

  private static Distribution d1(int lo, int hi, ReadGroupStats stats, boolean reverse) {
    final double offset = offsetLeft(stats, stats.fragmentMean() - stats.alignmentStartIgnored(), reverse);
    return DistributionUtils.distribution(lo, hi, stats.properRate(), offset, stats.fragmentStdDev(), reverse);
  }
  private static Distribution d2(int lo, int hi, ReadGroupStats stats, boolean reverse) {
    final double offset = offsetLeft(stats, stats.fragmentMean() - stats.meanLength() + stats.alignmentStartIgnored() + 1, reverse);
    return DistributionUtils.distribution(lo, hi, stats.properRate(), offset, stats.fragmentStdDev(), !reverse);
  }

  BreakpointBayesianSignal(final boolean debug) {
    super(debug);
  }

  BreakpointBayesianSignal() {
    super();
  }

  @Override
  Distribution leftArmProper(ReadGroupStats stats, boolean reverse) {
    final int lo = stats.lo();
    final int hi = stats.hi();
    final Distribution d1 = d1(lo, hi, stats, reverse);
    final double offset = offsetLeft(stats, stats.alignmentStartIgnored() + 1, reverse);
    final Distribution s = new DistributionStep(lo, hi, (int) offset, 0.0, stats.properRate(), reverse);
    //System.err.println("reverse=" + reverse + " a=" + stats.maxAlignment() + " offset=" + offset + " r=" + stats.meanLength() + " s[offset-1]=" + s.getValue((int) offset - 1) + " s[offset]=" + s.getValue((int) offset) + " s[offset+1]=" + s.getValue((int) offset + 1));
    final Distribution m = DistributionUtils.add(d1, s);
    return DistributionUtils.add(m, stats.properRandomRate());
  }
  @Override
  Distribution leftArmDiscordant(ReadGroupStats stats, boolean reverse) {
    final int lo = stats.lo();
    final int hi = stats.hi();
    final Distribution d2 = d2(lo, hi, stats, reverse);
    final double offset = offsetLeft(stats, stats.meanLength() - stats.alignmentStartIgnored(), reverse);
    final Distribution s = new DistributionStep(lo, hi, (int) offset, 1.0, 0.0, reverse);
    final Distribution m = DistributionUtils.multiply(d2, s);
    return DistributionUtils.add(m, stats.discordantRate());
  }
  @Override
  Distribution leftArmUnmated(ReadGroupStats stats, boolean reverse) {
    final int lo = stats.lo();
    final int hi = stats.hi();
    final Distribution d1 = d1(lo, hi, stats, reverse);
    final Distribution d2 = d2(lo, hi, stats, reverse);
    final Distribution a = DistributionUtils.add(d1, d2);
    final Distribution s = DistributionUtils.subtract(stats.properRate(), a);
    return DistributionUtils.add(s, stats.unmatedRate());
  }
  @Override
  Distribution rightArmProper(ReadGroupStats stats, boolean reverse) {
    return leftArmProper(stats, !reverse);
  }
  @Override
  Distribution rightArmDiscordant(ReadGroupStats stats, boolean reverse) {
    return leftArmDiscordant(stats, !reverse);
  }
  @Override
  Distribution rightArmUnmated(ReadGroupStats stats, boolean reverse) {
    return leftArmUnmated(stats, !reverse);
  }
}

