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

    assertFalse(r.completeGenomics());
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
