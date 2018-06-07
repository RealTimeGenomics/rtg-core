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

