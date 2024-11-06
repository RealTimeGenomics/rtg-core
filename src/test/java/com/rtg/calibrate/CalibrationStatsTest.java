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
