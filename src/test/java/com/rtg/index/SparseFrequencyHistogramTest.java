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

package com.rtg.index;

import junit.framework.TestCase;

/**
 */
public class SparseFrequencyHistogramTest extends TestCase {

  public void testBasic() {
    final SparseFrequencyHistogram hist = new SparseFrequencyHistogram();
    assertEquals(0, hist.length());
    try {
      hist.getFrequency(0);
      fail();
    } catch (ArrayIndexOutOfBoundsException e) {
      assertEquals("Array index out of range: 0", e.getMessage());
    }
    try {
      hist.getCount(-1);
      fail();
    } catch (ArrayIndexOutOfBoundsException e) {
      assertEquals("Array index out of range: -1", e.getMessage());
    }
    hist.add(2, 3);
    assertEquals(1, hist.length());
    assertEquals(2, hist.getFrequency(0));
    assertEquals(3, hist.getCount(0));
    try {
      hist.add(1, 1);
      fail();
    } catch (IllegalArgumentException e) {
      assertEquals("Frequencies must be added in ascending order.", e.getMessage());
    }
    hist.add(2, 4);
    assertEquals(1, hist.length());
    assertEquals(2, hist.getFrequency(0));
    assertEquals(7, hist.getCount(0));
  }

  public void testStupidPoints() {
    final SparseFrequencyHistogram hist = new SparseFrequencyHistogram();
    assertEquals(0, hist.arrayLengths());
    hist.add(Integer.MIN_VALUE, 1);
    hist.add(Integer.MIN_VALUE, 1);
  }

  public void testCalculateNewArrayLength() {
    assertEquals(1, SparseFrequencyHistogram.calculateNewArrayLength(0, 1));
    assertEquals(Integer.MAX_VALUE, SparseFrequencyHistogram.calculateNewArrayLength(Integer.MAX_VALUE / 3 * 2 + 2, Integer.MAX_VALUE / 3 * 2 + 3));
    assertEquals(15, SparseFrequencyHistogram.calculateNewArrayLength(10, 11));
  }

  public void testFromIndividualFrequencies() {
    final int[] freqs = {2, 3, 2, 5, 1, 1, 2, 1, 3, 4, 0, 0, 0};
    final SparseFrequencyHistogram hist = SparseFrequencyHistogram.fromIndividualFrequencies(freqs, 10);
    assertEquals(5, hist.length());
    assertEquals(1, hist.getFrequency(0));
    assertEquals(3, hist.getCount(0));
    assertEquals(2, hist.getFrequency(1));
    assertEquals(3, hist.getCount(1));
    assertEquals(3, hist.getFrequency(2));
    assertEquals(2, hist.getCount(2));
    assertEquals(4, hist.getFrequency(3));
    assertEquals(1, hist.getCount(3));
    assertEquals(5, hist.getFrequency(4));
    assertEquals(1, hist.getCount(4));
  }

  public void testMerge() {
    final SparseFrequencyHistogram histA = new SparseFrequencyHistogram();
    histA.add(1, 10);
    histA.add(3, 30);
    histA.add(4, 40);
    histA.add(5, 50);
    histA.add(7, 70);
    histA.add(8, 80);
    histA.add(9, 90);
    final SparseFrequencyHistogram histB = new SparseFrequencyHistogram();
    histB.add(2, 20);
    histB.add(4, 40);
    histB.add(6, 60);
    histB.add(8, 80);

    final SparseFrequencyHistogram histAB = SparseFrequencyHistogram.merge(histA, histB);
    final SparseFrequencyHistogram histBA = SparseFrequencyHistogram.merge(histB, histA);

    assertEquals(histAB.length(), histBA.length());
    for (int i = 0; i < histAB.length(); ++i) {
      assertEquals(histAB.getFrequency(i), histBA.getFrequency(i));
      assertEquals(histAB.getCount(i), histBA.getCount(i));
    }

    assertEquals(9, histAB.length());
    assertEquals(1, histAB.getFrequency(0));
    assertEquals(10, histAB.getCount(0));
    assertEquals(2, histAB.getFrequency(1));
    assertEquals(20, histAB.getCount(1));
    assertEquals(3, histAB.getFrequency(2));
    assertEquals(30, histAB.getCount(2));
    assertEquals(4, histAB.getFrequency(3));
    assertEquals(80, histAB.getCount(3));
    assertEquals(5, histAB.getFrequency(4));
    assertEquals(50, histAB.getCount(4));
    assertEquals(6, histAB.getFrequency(5));
    assertEquals(60, histAB.getCount(5));
    assertEquals(7, histAB.getFrequency(6));
    assertEquals(70, histAB.getCount(6));
    assertEquals(8, histAB.getFrequency(7));
    assertEquals(160, histAB.getCount(7));
    assertEquals(9, histAB.getFrequency(8));
    assertEquals(90, histAB.getCount(8));
  }
}
