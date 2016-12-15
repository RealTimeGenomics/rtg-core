/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.cnv.preprocess;

import java.util.Random;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class RegionDatasetTest extends TestCase {

  private RegionDataset makeDataSet(final int[] lengths, final double[] value) {
    final Random r = new Random();
    assertEquals(lengths.length, value.length);
    final RegionDataset ds = new RegionDataset(new String[] {"numeric-data"});
    ds.asNumeric(0);
    for (int k = 0; k < lengths.length; ++k) {
      final int start = r.nextInt(1000000);
      ds.add("chr", start, start + lengths[k], value[k]);
    }
    return ds;
  }

  public void testWeightedMedian() {
    assertEquals(10.0, makeDataSet(new int[] {1, 1, 1}, new double[] {1, 10, 100}).weightedMedian(0));
    assertEquals(100.0, makeDataSet(new int[] {1, 1, 3}, new double[] {1, 10, 100}).weightedMedian(0));
    assertEquals(100.0, makeDataSet(new int[] {1, 2, 3}, new double[] {1, 10, 100}).weightedMedian(0));
    assertEquals(10.0, makeDataSet(new int[] {1, 3, 3}, new double[] {1, 10, 100}).weightedMedian(0));
    assertEquals(100.0, makeDataSet(new int[] {1, 2, 2, 4}, new double[] {1, 10, 100, 1000}).weightedMedian(0));
  }

  public void testWeightedMean() {
    assertEquals(37.0, makeDataSet(new int[] {1, 1, 1}, new double[] {1, 10, 100}).weightedMean(0), 1e-4);
    assertEquals(62.2, makeDataSet(new int[] {1, 1, 3}, new double[] {1, 10, 100}).weightedMean(0), 1e-4);
    assertEquals(53.5, makeDataSet(new int[] {1, 2, 3}, new double[] {1, 10, 100}).weightedMean(0), 1e-4);
    assertEquals(469.0, makeDataSet(new int[] {1, 2, 2, 4}, new double[] {1, 10, 100, 1000}).weightedMean(0), 1e-4);
  }

  public void testFilter() {
    final RegionDataset d = makeDataSet(new int[] {1, 2, 2, 4}, new double[] {1, 10, 100, 1000});
    final NumericColumn dc = d.asNumeric(0);
    final RegionDataset d2 = d.filter(row -> dc.get(row) > 9 && dc.get(row) < 200);
    assertEquals(d.columns(), d2.columns());
    for (int i = 0; i < d.columns(); ++i) {
      assertEquals(d.columnName(i), d2.columnName(i));
    }
    assertEquals(2, d2.size());
  }
}
