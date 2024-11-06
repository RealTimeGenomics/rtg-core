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

import junit.framework.TestCase;


/**
 */
public class HeterozygousBayesianSignalTest extends TestCase {

  public void testCombine() {
    final Distribution d1 = new DistributionStep(-1, 2, 0, 4, 2);
    final Distribution d2 = new DistributionStep(-1, 2, 0, 5, 13);

    final HeterozygousBayesianSignal hbs = new HeterozygousBayesianSignal(null, null);

    final Distribution dc = hbs.combineDistributions(d1, d2);

    assertEquals(-1, dc.lo());
    assertEquals(2, dc.hi());

    assertEquals(4.5d, dc.get(-1));
    assertEquals(4.5d, dc.get(0));
    assertEquals(7.5d, dc.get(1));
  }

  public void testReturns() {
    final BayesianSignal bs1 = new NormalBayesianSignal(2);
    final BayesianSignal bs2 = new NormalBayesianSignal(5);

    final HeterozygousBayesianSignal hbs = new HeterozygousBayesianSignal(bs1, bs2);

    final ReadGroupStats rgs = new ReadGroupStats("blah", 2, 2, 2, 2, 2, 2, 2, 2, 2, 2);
    assertNotNull(hbs.leftArmDiscordant(rgs, false));
    assertNotNull(hbs.leftArmProper(rgs, false));
    assertNotNull(hbs.leftArmUnmated(rgs, false));
    assertNotNull(hbs.rightArmDiscordant(rgs, false));
    assertNotNull(hbs.rightArmProper(rgs, false));
    assertNotNull(hbs.rightArmUnmated(rgs, false));
  }
}
