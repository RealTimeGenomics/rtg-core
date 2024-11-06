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
