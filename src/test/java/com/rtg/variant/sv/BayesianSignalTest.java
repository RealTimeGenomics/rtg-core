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
