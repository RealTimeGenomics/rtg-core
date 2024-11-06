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
