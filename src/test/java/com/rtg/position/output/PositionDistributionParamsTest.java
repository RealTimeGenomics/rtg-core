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

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class PositionDistributionParamsTest extends TestCase {

  PositionDistributionParams getParams(final double sub, final double indel, final Integer maxGap) {
    return new PositionDistributionParams(sub, indel, maxGap, 0);
  }

  PositionDistributionParams getParams(final Integer maxGap) {
    return new PositionDistributionParams(maxGap);
  }

  public void test1() {

    final PositionDistributionParams a1 = getParams(0.1, 0.01, 1);
    a1.integrity();
    assertEquals(1, (int) a1.maxGap());
    assertEquals(0.1, a1.subsProb());
    assertEquals(0.01, a1.indelOpenProb());
    assertEquals(0.01, a1.indelExtendProb());

    final PositionDistributionParams a2 = getParams(0.1, 0.01, 1);
    a2.integrity();
    final PositionDistributionParams b = getParams(0.1, 0.01, 2);
    b.integrity();
    final PositionDistributionParams c = getParams(0.1, 0.09, 1);
    c.integrity();
    final PositionDistributionParams e = getParams(0.2, 0.01, 1);
    e.integrity();
    TestUtils.equalsHashTest(new PositionDistributionParams[][] {{a1, a2}, {b}, {c}, {e}});
    assertEquals("subs=0.2 indelOpen=0.01 indelExtend=0.01 maxGap=1", e.toString());
  }

  public void test2() {

    final PositionDistributionParams a1 = getParams(1);
    a1.integrity();
    assertEquals(1, (int) a1.maxGap());
    assertTrue(Double.isNaN(a1.subsProb()));
    assertTrue(Double.isNaN(a1.indelOpenProb()));
    assertTrue(Double.isNaN(a1.indelExtendProb()));

    final PositionDistributionParams a2 = getParams(1);
    a2.integrity();
    final PositionDistributionParams b = getParams(2);
    b.integrity();
    final PositionDistributionParams c = getParams(0.1, 0.09, 1);
    c.integrity();
    TestUtils.equalsHashTest(new PositionDistributionParams[][] {{a1, a2}, {b}, {c}});
    assertEquals("maxGap=1", a1.toString());
  }

}

