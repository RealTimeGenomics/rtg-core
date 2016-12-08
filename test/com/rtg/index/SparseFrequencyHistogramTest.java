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
