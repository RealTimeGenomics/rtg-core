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
