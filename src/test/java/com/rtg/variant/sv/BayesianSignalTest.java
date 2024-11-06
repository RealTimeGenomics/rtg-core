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

import com.rtg.sam.ReadGroupUtils;

import junit.framework.TestCase;

/**
 */
public class BayesianSignalTest extends TestCase {

  public void testCompactConstant() {
    final DistributionConstant d1 = new DistributionConstant(-50, 50, 33.3);
    final DistributionConstant d2 = (DistributionConstant) BayesianSignal.compactDistribution(d1, "");

    assertEquals(d1.lo(), d2.lo());
    assertEquals(d1.hi(), d2.hi());
    assertEquals(d1.getConstant(), d2.getConstant());

  }

  public void testCompactStep() {
    final DistributionStep d1 = new DistributionStep(-50, 50, 3, 33.3, 45.5);
    final DistributionStep d2 = (DistributionStep) BayesianSignal.compactDistribution(d1, "");

    assertEquals(d1.lo(), d2.lo());
    assertEquals(d1.hi(), d2.hi());
    assertEquals(d1.getOffset(), d2.getOffset());
    assertEquals(d1.getRate1(), d2.getRate1());
    assertEquals(d1.getRate2(), d2.getRate2());
  }

  public void testOffset() {
    final ReadGroupStats stats = new ReadGroupStats(ReadGroupUtils.UNKNOWN_RG, 20, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0);
    assertEquals(0.0, BayesianSignal.offsetLeft(stats, 0.0, false), 0.000001);
    assertEquals(-21.0, BayesianSignal.offsetLeft(stats, 0.0, true));

    assertEquals(-5.0, BayesianSignal.offsetLeft(stats, 5.0, false));
    assertEquals(-16.0, BayesianSignal.offsetLeft(stats, 5.0, true));
  }
}
