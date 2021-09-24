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
package com.rtg.ngs;

import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class NgsFilterParamsTest extends TestCase {

  NgsFilterParams getParams(final int topN, final int errorLimit, final boolean exclude, final boolean useIds, final boolean zip) {
    return NgsFilterParams.builder().outputFilter(OutputFilter.NONE).zip(zip).topN(topN).exclude(exclude).useids(useIds).errorLimit(errorLimit).unmatedMaxMismatches(IntegerOrPercentage.valueOf(4)).matedMaxMismatches(IntegerOrPercentage.valueOf(7)).create();
  }

  public void testDefaultAS() {
    final NgsFilterParams a1 = NgsFilterParams.builder().topN(1).zip(false).create();
    assertNotNull(a1.matedMaxMismatches());
  }

  public void testEquals() {
    final NgsFilterParams a1 = getParams(10, 5, false, false, false);
    final NgsFilterParams a2 = getParams(10, 5, false, false, false);
    final NgsFilterParams b = getParams(10, 5, false, false, true);
    final NgsFilterParams d = getParams(10, 3, false, false, false);
    final NgsFilterParams e = getParams(10, 5, false, true, false);
    final NgsFilterParams f = getParams(11, 5, false, false, false);
    final NgsFilterParams g = getParams(10, 5, true, false, false);
    final NgsFilterParams h = getParams(10, 5, true, true, false);
    TestUtils.equalsHashTest(new NgsFilterParams[][] {{a1, a2}, {b}, {d}, {e}, {f}, {g}, {h}});
  }

  public void test0() {
    final NgsFilterParams sp = getParams(10, 5, false, false, false);
    sp.integrity();
    assertEquals(5, sp.errorLimit());
    assertEquals(OutputFilter.NONE, sp.outputFilter());
    assertEquals(10, sp.topN());
    assertEquals(false, sp.exclude());
    assertEquals(10.0, sp.maxEScore());
    assertEquals(0.0, sp.minBitScore());
    assertEquals(50, sp.minIdentity());
    assertEquals(50, sp.preFilterMinScore());
    assertEquals(50, sp.preFilterMinOverlap());
    assertEquals(-3, sp.preFilterAlgorithm());
    assertEquals(""
        + "NgsFilterParams filter=NONE topN=10 maxTopResults=5 error limit=5 mated max mismatches=7 unmated max mismatches=4 exclude=" + Boolean.FALSE.toString()
        + " use-ids=" + Boolean.FALSE.toString()
        + " zip=" + Boolean.FALSE.toString()
        , sp.toString()
    );
  }

  public void testSetPreFilter() {
    final NgsFilterParams fp = NgsFilterParams.builder().outputFilter(OutputFilter.NONE)
    .topN(10).errorLimit(5).unmatedMaxMismatches(IntegerOrPercentage.valueOf(4))
    //.matedMaxMismatches(IntegerOrPercentage.valueOf(7))
    .preFilterAlgorithm(2).preFilterMinScore(0).preFilterMinOverlap(100)
    .create();
    assertEquals(2, fp.preFilterAlgorithm());
    assertEquals(0, fp.preFilterMinScore());
    assertEquals(100, fp.preFilterMinOverlap());
  }
}
