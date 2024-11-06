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
package com.rtg.position.output;

import junit.framework.TestCase;

/**
 * See the spreadsheet for the test case below.
 */
public class GappedDistributionTest extends TestCase {

  public void test() {
    final GappedDistribution gd = new GappedDistribution(8, 8, GappedDistribution.distrParams(16));
    gd.integrity();
    final String s = gd.toString();
    //System.err.println(s);
    assertTrue(s.startsWith("[ 24]= "));
    assertTrue(s.contains("     i "));
    assertTrue(s.contains("     d "));
    assertTrue(s.contains("     t "));
    assertTrue(s.contains(" * 0.0098")); //total at end of first line
    assertTrue(s.contains(" * 1.0045")); //total at end of last line
    assertTrue(s.contains("     t 1.0000= 0.0045")); //total first entry line 0
    assertTrue(s.contains("     t 0.0045  0.9910= 0.0089")); //total first entry line 1
    assertTrue(s.contains("0.0012  0.0045  0.0023= 0.0006  0.0011  0.0001")); //values at maximum for row 16
    assertTrue(s.contains("0.0005  0.0296  0.0083= 0.0335  0.0008")); //values at maximum for row 8
    gd.probabilities(); //make sure doesn't fall over - really tested in GapProbabilitiesTest
  }

  public void testProb() {
    assertTrue(GappedDistribution.prob(Double.MIN_VALUE));
    assertTrue(GappedDistribution.prob(1.0));
    assertTrue(GappedDistribution.prob(0.5));

    failProb(0.0);
    failProb(1.00001);
    failProb(Double.NaN);
    failProb(Double.POSITIVE_INFINITY);
    failProb(Double.NEGATIVE_INFINITY);
  }

  private void failProb(final double x) {
    try {
      GappedDistribution.prob(x);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("" + x, e.getMessage());
    }
  }

  public void testStepOffset() {
    assertEquals(0, GappedDistribution.stepOffset(1, 1));
    assertEquals(0, GappedDistribution.stepOffset(3, 9));
    assertEquals(0, GappedDistribution.stepOffset(3, 3));
    assertEquals(1, GappedDistribution.stepOffset(3, 5));
    assertEquals(4, GappedDistribution.stepOffset(8, 12));
  }
}
