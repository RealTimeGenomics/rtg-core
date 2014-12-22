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

package com.rtg.util;

import junit.framework.TestCase;

/**
 */
public class HistogramTest extends TestCase {

  public void test() {
    final Histogram hist = new Histogram();
    assertEquals(0, hist.getLength());
    assertEquals("", hist.toString());
    hist.increment(3);
    assertEquals(4, hist.getLength());
    assertEquals("0\t0\t0\t1", hist.toString());
    for (int i = 0; i < 3; i++) {
      assertEquals(0, hist.getValue(i));
    }
    assertEquals(1, hist.getValue(3));
    hist.increment(2, 10);
    assertEquals(4, hist.getLength());
    assertEquals(10, hist.getValue(2));
    hist.increment(0);
    assertEquals("1\t0\t10\t1", hist.toString());
    hist.increment(9, 9);
    assertEquals(10, hist.getLength());
    assertEquals("1\t0\t10\t1\t0\t0\t0\t0\t0\t9", hist.toString());
    final Histogram hist2 = new Histogram();
    hist2.addHistogram(hist.toString());
    assertEquals(10, hist2.getLength());
    assertEquals("1\t0\t10\t1\t0\t0\t0\t0\t0\t9", hist2.toString());
    final double[] dist = hist.toDistribution();
    assertEquals(10, dist.length);
    assertEquals(0.0476, dist[0], 0.00005);
    assertEquals(0.0, dist[1], 0.00005);
    assertEquals(0.4762, dist[2], 0.00005);
    assertEquals(0.0476, dist[3], 0.00005);
    assertEquals(0.0, dist[4], 0.00005);
    assertEquals(0.0, dist[5], 0.00005);
    assertEquals(0.0, dist[6], 0.00005);
    assertEquals(0.0, dist[7], 0.00005);
    assertEquals(0.0, dist[8], 0.00005);
    assertEquals(0.4286, dist[9], 0.00005);
  }

  public void testEmptyParse() {
    final Histogram hist = new Histogram();
    hist.addHistogram("");
    assertEquals(0, hist.getLength());
    assertEquals("", hist.toString());
  }
}
