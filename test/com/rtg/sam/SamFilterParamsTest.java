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
package com.rtg.sam;

import com.rtg.util.IntegerOrPercentage;

import junit.framework.TestCase;

/**
 * Test class.
 */
public class SamFilterParamsTest extends TestCase {

  public void testDefaults() {
    // Default behaviour should be to not filter any records.
    final SamFilterParams.SamFilterParamsBuilder builder = SamFilterParams.builder();
    assertNotNull(builder);
    final SamFilterParams params = builder.create();
    assertNotNull(params);
    assertEquals(-1, params.maxAlignmentCount());
    assertEquals(null, params.maxMatedAlignmentScore());
    assertEquals(null, params.maxUnmatedAlignmentScore());
    assertEquals(-1, params.minMapQ());
    assertFalse(params.findAndRemoveDuplicates());
    assertEquals(0, params.requireSetFlags());
    assertEquals(0, params.requireUnsetFlags());
    assertFalse(params.excludeUnmated());
  }

  public void testActual() {
    final SamFilterParams p = SamFilterParams.builder()
      .minMapQ(10)
      .maxAlignmentCount(42)
      .maxMatedAlignmentScore(new IntegerOrPercentage(43))
      .maxUnmatedAlignmentScore(new IntegerOrPercentage(44))
      .excludeMated(true)
      .excludeUnmated(true)
      .excludeUnmapped(false)
      .excludeDuplicates(true)
      .create();
    assertEquals("SamFilterParams minMapQ=10 maxAlignmentCount=42 maxMatedAlignmentScore=43 maxUnmatedAlignmentScore=44 excludeUnmated=" + true + " excludeUnplaced=" + false + " requireSetFlags=0 requireUnsetFlags=1026 regionTemplate=" + null, p.toString());
    assertEquals(42, p.maxAlignmentCount());
    assertEquals(new IntegerOrPercentage(43), p.maxMatedAlignmentScore());
    assertEquals(new IntegerOrPercentage(44), p.maxUnmatedAlignmentScore());
    assertEquals(1026, p.requireUnsetFlags());
    assertEquals(0, p.requireSetFlags());
    assertTrue(p.excludeUnmated());
  }

  public void testActual2() {
    final SamFilterParams p = SamFilterParams.builder()
      .maxUnmatedAlignmentScore(new IntegerOrPercentage(44))
      .excludeUnmated(true)
      .excludeUnmapped(true)
      .excludeUnplaced(true)
      .requireSetFlags(256)
      .create();
    assertEquals("SamFilterParams minMapQ=-1 maxAlignmentCount=-1 maxMatedAlignmentScore=null maxUnmatedAlignmentScore=44 excludeUnmated=" + true + " excludeUnplaced=" + true + " requireSetFlags=256 requireUnsetFlags=4 regionTemplate=" + null, p.toString());
    assertEquals(-1, p.maxAlignmentCount());
    assertEquals(null, p.maxMatedAlignmentScore());
    assertEquals(new IntegerOrPercentage(44), p.maxUnmatedAlignmentScore());
    assertEquals(256, p.requireSetFlags());
    assertTrue(p.excludeUnmated());
  }

}
