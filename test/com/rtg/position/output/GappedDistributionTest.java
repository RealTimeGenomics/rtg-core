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

import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * See the spreadsheet for the test case below.
 */
public class GappedDistributionTest extends TestCase {

  public static TestSuite suite() {
    final TestSuite suite = new TestSuite();
    suite.addTestSuite(GapProbabilitiesScorerTest.class);
    suite.addTestSuite(GappedDistributionTest.class);
    return suite;
  }


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
