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
package com.rtg.variant.util;

import java.io.IOException;

import junit.framework.TestCase;

/**
 */
public class GenotypeTest extends TestCase {

  public void testDip() throws IOException {
    Genotype gt = new Genotype(new int[] {1, 0});
    assertEquals(2, gt.length());
    assertEquals(0, gt.get(0));
    assertEquals(1, gt.get(1));
    assertEquals("0/1", gt.toString());

    assertTrue(gt.contains(0));
    assertTrue(gt.contains(1));
    assertFalse(gt.contains(2));
  }

  public void testHap() throws IOException {
    Genotype gt = new Genotype(new int[]{2});
    assertEquals(1, gt.length());
    assertEquals(2, gt.get(0));
    assertEquals("2", gt.toString());

    assertFalse(gt.contains(0));
    assertFalse(gt.contains(1));
    assertTrue(gt.contains(2));
  }
}
