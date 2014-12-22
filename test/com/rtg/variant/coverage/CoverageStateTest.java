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
package com.rtg.variant.coverage;

import com.rtg.util.MathUtils;

import junit.framework.TestCase;

/**
 */
public class CoverageStateTest extends TestCase {

  public void test() {
    final CoverageState cs = new CoverageState("t1", new byte[] {1, 2, 3}, true);
    assertEquals("t1", cs.getName());
    assertEquals(2, cs.getTemplate(1));
    assertEquals(3, cs.getTemplateLength());
    assertFalse(cs.isSeen());
    assertTrue(cs.isAdditional());
    cs.setSeen();
    assertTrue(cs.isSeen());
    assertEquals(0.0, cs.getCoverage(1), 1e-4);
    cs.incrementCoverage(1, 1);
    assertEquals(1.0, cs.getCoverage(1), 1e-4);
    cs.incrementCoverage(1, 0.5);
    assertEquals(1.5, cs.getCoverage(1), 1e-4);
    cs.incrementCoverage(1, 50000);
    assertEquals(50001.5, cs.getCoverage(1), 1e-4);
    cs.incrementCoverage(1, 1.0 / 3.0);
    assertEquals(50001.5 + 1.0 / 3.0, cs.getCoverage(1), 1e-4);
    assertEquals(0, cs.getIH1(1));
    assertEquals(0, cs.getIHgt1(1));
    cs.incrementIH(1, 1);
    assertEquals(1, cs.getIH1(1));
    assertEquals(0, cs.getIHgt1(1));
    cs.incrementIH(1, 2);
    assertEquals(1, cs.getIH1(1));
    assertEquals(1, cs.getIHgt1(1));
    for (int k = 0; k < 10; k++) {
      cs.incrementCoverage(2, 1);
      cs.incrementCoverage(2, 0.5);
      cs.incrementCoverage(2, 0.5);
      cs.incrementCoverage(2, 1.0 / 3.0);
      cs.incrementCoverage(2, 1.0 / 3.0);
      cs.incrementCoverage(2, 1.0 / 3.0);
      cs.incrementCoverage(2, 0.25);
      cs.incrementCoverage(2, 0.25);
      cs.incrementCoverage(2, 0.25);
      cs.incrementCoverage(2, 0.25);
      cs.incrementCoverage(2, 0.2);
      cs.incrementCoverage(2, 0.2);
      cs.incrementCoverage(2, 0.2);
      cs.incrementCoverage(2, 0.2);
      cs.incrementCoverage(2, 0.2);
    }
    assertEquals(50, (int) MathUtils.round(cs.getCoverage(2)));
  }

  public void test2() {
    final CoverageState cs = new CoverageState("t1", new byte[] {1, 2, 3}, false);
    assertEquals("t1", cs.getName());
    assertFalse(cs.isAdditional());
    double t = 0;
    cs.incrementCoverage(2, 0.5);
    t += 0.5;
    cs.incrementCoverage(2, 1.0 / 3.0);
    t += 1.0 / 3.0;
    cs.incrementCoverage(2, 0.5);
    t += 0.5;
    cs.incrementCoverage(2, 1.0 / 3.0);
    t += 1.0 / 3.0;
    cs.incrementCoverage(2, 0.5);
    t += 0.5;
    cs.incrementCoverage(2, 1.0 / 3.0);
    t += 1.0 / 3.0;
    cs.incrementCoverage(2, 0.5);
    t += 0.5;
    cs.incrementCoverage(2, 1.0 / 3.0);
    t += 1.0 / 3.0;
    cs.incrementCoverage(2, 0.5);
    t += 0.5;
    cs.incrementCoverage(2, 1.0 / 3.0);
    t += 1.0 / 3.0;
    cs.incrementCoverage(2, 0.5);
    t += 0.5;
    cs.incrementCoverage(2, 1.0 / 3.0);
    t += 1.0 / 3.0;
    cs.incrementCoverage(2, 0.5);
    t += 0.5;
    cs.incrementCoverage(2, 1.0 / 3.0);
    t += 1.0 / 3.0;
    cs.incrementCoverage(2, 0.5);
    t += 0.5;
    cs.incrementCoverage(2, 1.0 / 3.0);
    t += 1.0 / 3.0;
    cs.incrementCoverage(2, 0.5);
    t += 0.5;
    cs.incrementCoverage(2, 1.0 / 3.0);
    t += 1.0 / 3.0;
    // 9 * (1/2) + 9 * (1/3) = 7.5, or 8 when rounded
    assertEquals(8, MathUtils.round(cs.getCoverage(2)));
    // Note this is a case that double arithmetic gets wrong ...
    assertEquals(7, MathUtils.round(t));
  }
}
