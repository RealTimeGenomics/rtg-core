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

package com.rtg.calibrate;

import junit.framework.TestCase;

/**
 */
public class CalibrationStatsTest extends TestCase {

  public void test() {
    final Covariate[] covs = {new CovariateReadGroup(), new CovariateSequence()};
    covs[0].parse("Foo");
    covs[1].parse("Bar");
    final CalibrationStats stats = new CalibrationStats(new int[] {0, 0});
    assertEquals("Foo\tBar", stats.getName(covs));
    assertEquals("Foo\tBar\t0\t0\t0\t0", stats.outputString(covs));
    stats.seenEquals(3);
    stats.seenEquals(5);
    stats.seenMnp(7);
    stats.seenMnp(11);
    stats.seenInsert(13);
    stats.seenDelete(17);
    assertEquals("Foo\tBar\t8\t18\t13\t17", stats.outputString(covs));
    assertEquals(56, stats.getTotalLength());
    assertEquals(8, stats.getEqual());
    assertEquals(18, stats.getDifferent());
    assertEquals(13, stats.getInserted());
    assertEquals(17, stats.getDeleted());
  }

  public void testAccumulate() {
    final Covariate[] covsA = {new CovariateReadGroup(), new CovariateSequence()};
    covsA[0].parse("Foo");
    covsA[1].parse("Bar");
    final CalibrationStats stats1 = new CalibrationStats(new int[] {0, 0});
    stats1.seenEquals(3);
    stats1.seenMnp(7);
    stats1.seenDelete(1);
    final Covariate[] covsB = {new CovariateReadGroup()};
    covsB[0].parse("Etc");
    final CalibrationStats stats2 = new CalibrationStats(new int[] {0});
    stats2.seenEquals(11);
    stats2.seenMnp(12);
    stats2.seenDelete(100);
    stats2.seenInsert(200);
    stats1.accumulate(stats2);
    assertEquals("Foo\tBar\t14\t19\t200\t101", stats1.outputString(covsA));
    assertEquals("Etc\t11\t12\t200\t100", stats2.outputString(covsB));
  }
}
