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
package com.rtg.index.similarity;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;



/**
 */
public class SimilarityMatrixTest extends TestCase {

  /**
   * Test method for {@link com.rtg.index.similarity.SimilarityMatrix}.
   */
  public final void test() {
    final SimilarityMatrix sim = new SimilarityMatrix(10);
    assertEquals(10, sim.length());
    sim.globalIntegrity();
    for (int i = 0; i < 10; ++i) {
      for (int j = 0; j < 10; ++j) {
        assertEquals(0, sim.get(i, j), 0);
      }
    }

    sim.increment(0, 0);

    sim.increment(2, 5);
    sim.increment(5, 2);

    sim.increment(9, 9);

    sim.increment(6, 2, 3);
    sim.increment(2, 6, 5);

    sim.globalIntegrity();

    assertEquals(1, sim.get(0, 0), 1.0E-8);
    assertEquals(1, sim.get(9, 9), 1.0E-8);

    assertEquals(2, sim.get(2, 5), 1.0E-8);
    assertEquals(2, sim.get(5, 2), 1.0E-8);

    assertEquals(8, sim.get(2, 6), 1.0E-8);
    assertEquals(8, sim.get(6, 2), 1.0E-8);

    for (int i = 0; i < 1001; ++i) {
      sim.increment(3, 4);
    }
    assertEquals(1001, sim.get(4, 3), 1.0E-8);
    final String expected = ""
      + "SimilarityMatrix 10" + StringUtils.LS
      + "[0]\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[1]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[2]\t0\t0\t0\t0\t0\t2\t8\t0\t0\t0" + StringUtils.LS
      + "[3]\t0\t0\t0\t0\t1001\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[4]\t0\t0\t0\t1001\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[5]\t0\t0\t2\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[6]\t0\t0\t8\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[7]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[8]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[9]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t1" + StringUtils.LS
      + StringUtils.LS
      ;
    assertEquals(expected, sim.toString());
  }

  public void testSet() {
    final SimilarityMatrix sim = new SimilarityMatrix(10);
    sim.globalIntegrity();

    assertEquals(0, sim.get(3, 4), 1.0E-8);
    assertEquals(0, sim.get(4, 3), 1.0E-8);

    sim.set(3, 4, 5);
    assertEquals(5, sim.get(3, 4), 1.0E-8);
    assertEquals(5, sim.get(4, 3), 1.0E-8);

    sim.set(4, 3, 0);
    assertEquals(0, sim.get(3, 4), 1.0E-8);
    assertEquals(0, sim.get(4, 3), 1.0E-8);

    try {
      sim.set(4, 3, -1);
      fail();
    } catch (final RuntimeException e) {
      //ok
    }

  }

  /**
   * Test method for {@link com.rtg.index.similarity.SimilarityMatrix}.
   */
  public final void testBigBad() {
    try {
      new SimilarityMatrix(Integer.MAX_VALUE + 1L);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("2147483648", e.getMessage());
    }
  }
}
