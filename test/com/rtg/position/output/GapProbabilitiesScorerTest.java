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

import static com.rtg.position.output.GapProbabilitiesScorer.PROB_FORMAT;

import junit.framework.TestCase;

/**
 */
public class GapProbabilitiesScorerTest extends TestCase {

  public void test() {
    //See GappedDistribution.xls for a spreadsheet which computes this case
    final GappedDistribution gd = new GappedDistribution(8, 8, GappedDistribution.distrParams(16));
    gd.integrity();
    final GapProbabilitiesScorer gp = (GapProbabilitiesScorer) gd.probabilities();
    gp.globalIntegrity();
    //System.err.println(gp.toString());
    assertEquals(-1, gp.minDelta());
    assertEquals(+1, gp.maxDelta());

    assertEquals(0.0, gp.score(0, -1, 0, -1));
    assertEquals("-5.4037", PROB_FORMAT.format(gp.score(0, 0, 0, 1)));
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 0, 0, 2));

    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 8, 0, 6));
    assertEquals("-3.3870", PROB_FORMAT.format(gp.score(0, 8, 0, 7)));
    assertEquals("-4.7621", PROB_FORMAT.format(gp.score(0, 8, 0, 8)));
    assertEquals("-3.3948", PROB_FORMAT.format(gp.score(0, 8, 0, 9)));
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 8, 0, 10));

    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 16, 0, 14));
    assertEquals("-5.3992", PROB_FORMAT.format(gp.score(0, 16, 0, 15)));
    assertEquals("-6.0065", PROB_FORMAT.format(gp.score(0, 16, 0, 16)));
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 16, 0, 17));

    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 24, 0, 24));
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 0, 0, 17));

    assertEquals(" 0.0000", PROB_FORMAT.format(gp.scoreMax(0, 0, -1, -1)));
    assertEquals("-0.0090", PROB_FORMAT.format(gp.scoreMax(0, 1, -1, -1)));
    assertEquals("-3.3870", PROB_FORMAT.format(gp.scoreMax(0, 8, -1, -1)));
    assertEquals("-3.3948", PROB_FORMAT.format(gp.scoreMax(0, 9, -1, -1)));
    assertEquals("-5.3992", PROB_FORMAT.format(gp.scoreMax(0, 16, -1, -1)));
    assertEquals("-5.4044", PROB_FORMAT.format(gp.scoreMax(0, 17, -1, -1)));
    assertEquals("-5.4200", PROB_FORMAT.format(gp.scoreMax(0, 20, -1, -1)));
    assertEquals(Double.NEGATIVE_INFINITY, gp.scoreMax(0, 25, -1, -1));

    assertEquals(0.0, gp.scoreMax(0, 0, -1, -1));
    assertEquals("-5.3992", PROB_FORMAT.format(gp.scoreMax(0, 16, -1, -1)));
    assertEquals(Double.NEGATIVE_INFINITY, gp.scoreMax(0, 25, -1, -1));
  }

  /** The word size is not a multiple of the stepsize so there is an offset when extracting the probabilities. */
  public void testOffset() {
    final GappedDistribution gd = new GappedDistribution(8, 12, GappedDistribution.distrParams(16));
    gd.integrity();
    final GapProbabilitiesScorer gp = (GapProbabilitiesScorer) gd.probabilities();
    gp.globalIntegrity();
    //System.err.println(gp.toString());
    assertEquals(-1, gp.minDelta());
    assertEquals(+1, gp.maxDelta());
    //System.err.println(gp.gapThreshold());

    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 4, 0, 2));
    assertEquals("-4.0444", PROB_FORMAT.format(gp.score(0, 4, 0, 3)));
    assertEquals("-0.0359", PROB_FORMAT.format(gp.score(0, 4, 0, 4)));
    assertEquals("-3.8303", PROB_FORMAT.format(gp.score(0, 4, 0, 5)));
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 4, 0, 6));

    assertEquals("-4.0444", PROB_FORMAT.format(gp.score(0, 4, 0, 3)));
    assertEquals("-0.0359", PROB_FORMAT.format(gp.score(0, 4, 0, 4)));
    assertEquals("-3.8303", PROB_FORMAT.format(gp.score(0, 4, 0, 5)));

    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 12, 0, 10));
    assertEquals("-3.4179", PROB_FORMAT.format(gp.score(0, 12, 0, 11)));
    assertEquals("-4.6649", PROB_FORMAT.format(gp.score(0, 12, 0, 12)));
    assertEquals("-3.4255", PROB_FORMAT.format(gp.score(0, 12, 0, 13)));
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 12, 0, 14));

    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 20, 0, 18));
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 20, 0, 19));
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 20, 0, 20));
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(0, 4, 0, 17));

    final String str = gp.toString();
    //System.err.println(str);
    assertTrue(str.startsWith("       [  0]  [  1]"));
    assertTrue(str.contains("[  0]  0.000 -5.404 - "));
    assertTrue(str.contains("[  3] -      -      -4.323 -0.027 -4.044 - "));
    assertTrue(str.contains("[  4] -      -      -      -4.044 -0.036 -3.830 - "));
    assertTrue(str.contains("[  5] -      -      -      -      -3.830 -0.045 -3.657 - "));
    assertTrue(str.contains("[ 11] -      -      -      -      -      -      -      -      -      -      -3.410 -4.688 -3.418 - "));
    assertTrue(str.contains("[ 12] -      -      -      -      -      -      -      -      -      -      -      -3.418 -4.665 -3.425 - "));
  }

  public void test1() {
    check(7, 19);
    check(1, 1);
    check(3, 5);
    check(3, 3);
  }

  private void check(final int st, final int gap) {
    check(st, st, gap);
  }

  public void test2() {
    check(3, 5, 19);
    check(8, 16, 16);
  }

  /**
   * Check when word size and stepsize differ
   */
  private void check(final int st, final int ws, final int gap) {
    final GappedDistribution gd = new GappedDistribution(st, ws, GappedDistribution.distrParams(gap));
    gd.integrity();
    final GapProbabilitiesScorer gp = (GapProbabilitiesScorer) gd.probabilities();
    gp.globalIntegrity();
    //System.err.println(gp.toString());
    if (ws % st == 0) {
      assertEquals(0.0, gp.score(0, -1, 0, -1));
    }
    final int offset = GappedDistribution.stepOffset(st, ws);
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(((gap - offset) / st + 1) * st + offset, -1, 0, -1));
    assertEquals(Double.NEGATIVE_INFINITY, gp.score(offset, -1, gap + 1, -1));
  }

  /**
   * Use a weird underlyng distrbution to force the minDelta to be negative and for that to occur on the first row.
   * Also maxDelta on first row.
   */
  public void testDelta() {
    final GapProbabilitiesScorer gp = new GapProbabilitiesScorer(0, 0.1, new double[][]{{0.0, 0.1, 0.1}}, new double[] {0.1});
    gp.globalIntegrity();
    assertEquals(-2, gp.minDelta());
    assertEquals(-1, gp.maxDelta());
  }

  /**
   * Use a weird underlying distrbution to force the maxDelta to be negative.
   */
  public void testMaxDelta() {
    final GapProbabilitiesScorer gp = new GapProbabilitiesScorer(0, 0.1, new double[][]{{0.0, 0.0, 0.0}, {0.0, 0.1, 0.1}}, new double[] {0.0, 0.1});
    gp.globalIntegrity();
    assertEquals(-1, gp.minDelta());
    assertEquals(0, gp.maxDelta());
  }

  /**
   * Use an all zero underlying distrbution to force the maxDelta and minDelta to be very large/small.
   */
  public void testEmptyDelta() {
    final GapProbabilitiesScorer gp = new GapProbabilitiesScorer(0, 0.1, new double[][]{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}, new double[] {0.0, 0.0});
    gp.globalIntegrity();
    assertEquals(Integer.MAX_VALUE, gp.minDelta());
    assertEquals(Integer.MIN_VALUE, gp.maxDelta());
  }

}
