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

package com.rtg.variant.realign;


import junit.framework.TestCase;

/**
 */
public class RealignParamsGenomeTest extends TestCase {
  public void test() {
    final RealignParamsGenome r = RealignParamsGenome.SINGLETON;
    /* Pointless checking of magic numbers
    assertEquals(-7.505592279737757, r.insertOpenLn(), 1e-8);
    assertEquals(-7.505592279737757, r.deleteOpenLn(), 1e-8);
    assertEquals(-0.6931471805599453, r.insertExtendLn(), 1e-8);
    assertEquals(-0.6931471805599453, r.deleteExtendLn(), 1e-8);
    assertEquals(0.00124, r.misMatch(), 1e-8);
    */
    assertEquals(Math.log(r.misMatch()), r.misMatchLn(), 1e-8);
    assertEquals(Math.log(1.0 - r.misMatch()), r.matchLn(), 1e-8);

    assertNull(r.machineType());
    try {
      r.gapStart(0);
      fail();
    } catch (final UnsupportedOperationException e) {
    }
    try {
      r.gapEnd(0);
      fail();
    } catch (final UnsupportedOperationException e) {
    }
    try {
      r.gapFreqLn(0, 0);
      fail();
    } catch (final UnsupportedOperationException e) {
    }
    try {
      r.gapDistributionPoss(null);
      fail();
    } catch (final UnsupportedOperationException e) {
    }

  }

}
