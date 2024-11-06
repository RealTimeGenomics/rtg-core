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
public class SignalDistributionLnTest extends TestCase {

  public void test() {
    final SamCounts sa = new SamArray(5);
    sa.increment(1);
    sa.increment(2);
    sa.increment(2);
    sa.increment(3);
    sa.increment(3);
    sa.increment(3);
    sa.increment(4);
    sa.increment(4);
    sa.increment(4);
    sa.increment(4);
    final SignalDistributionLn sig = new SignalDistributionLn(sa, new DistributionConstant(-2, 2, 2.0), "blah");
    sig.globalIntegrity();

    assertEquals(6.3069, sig.value(0), 0.0001);
    assertEquals(4.3069, sig.value(1), 0.0001);
    assertEquals(2.5232, sig.value(2), 0.0001);
    assertEquals(1.2958, sig.value(3), 0.0001);
    assertEquals(2.9890, sig.value(4), 0.0001);
    assertEquals("blah", sig.columnLabel());
  }

  //ensure that minimum signal occurs at the right place and is zero (exactly mimics distribution)
  public void testMin() {
    //distribution is a step 100 on left and 1 on right - actual signal has this distribution centered at 10
    final SamArray sa = new SamArray(20);
    for (int i = 0; i < 10; ++i) {
      sa.increment(i, 100);
    }
    for (int i = 10; i < 20; ++i) {
      sa.increment(i, 1);
    }
    final Distribution distr = new DistributionStep(-5, 5, -1, 100.0, 1.0);
    //System.err.println(distr.toString());
    //System.err.println(distr.dump());
    final int expMini = 10;
    final double expMin = 0.0;

    DistributionTestUtils.checkMin(sa, distr, expMini, expMin);
  }

  public void testMinErf() {
    //distribution is an erf distribution on left and 1 on right - actual signal has this distribution centered at 10
    final SamArray sa = new SamArray(20);
    for (int i = 0; i <= 4; ++i) {
      sa.increment(i, 101);
    }
    sa.increment(6, 100);
    sa.increment(6, 99);
    sa.increment(7, 94);
    sa.increment(8, 85);
    sa.increment(9, 70);
    sa.increment(10, 51);
    sa.increment(11, 32);
    sa.increment(12, 17);
    sa.increment(13, 8);
    sa.increment(14, 3);
    sa.increment(15, 2);
    for (int i = 16; i < 20; ++i) {
      sa.increment(i, 1);
    }
    final Distribution distr = DistributionUtils.add(DistributionUtils.distribution1(-3, 3, 100.0, 0.0, 2.0), 1.0);

    //System.err.println(distr.toString());
    //System.err.println(distr.dump());
    DistributionTestUtils.checkMin(sa, distr, 10, 0.00167);
  }

}
