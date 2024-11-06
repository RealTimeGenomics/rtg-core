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

package com.rtg.variant.cnv.preprocess;

import java.util.Arrays;
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

  public void testRestrictedColumns() {
    final String[] desiredColumns = {"second-data"};
    final RegionDataset d = new RegionDataset(desiredColumns);
    d.setColumnIndex(Arrays.asList("numeric-data", "second-data"));
    d.add("chr", 1, 2, 10, 110);
    d.add("chr", 3, 4, "12", "130");

    final NumericColumn dc = d.asNumeric(0);
    assertEquals(240.0, dc.sum());
    assertEquals(Arrays.asList(desiredColumns), d.getColumnNames());
  }

  public void testFilterIntColumn() {
    final RegionDataset d = new RegionDataset(new String[0]);
    d.addColumn(new IntColumn("mycol"));
    assertEquals(1, d.getColumns().size());
    final RegionDataset f = d.filter(row -> false);
    assertEquals(1, f.getColumns().size());
    assertTrue(f.getColumns().get(0) instanceof IntColumn);
  }
}
