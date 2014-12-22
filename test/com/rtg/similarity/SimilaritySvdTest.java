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

package com.rtg.similarity;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.GregorianCalendar;
import java.util.Random;

import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 * Tests for SimilaritySvd class.
 */
public class SimilaritySvdTest extends TestCase {
  private NanoRegression mNano;

  @Override
  public void setUp() throws Exception {
    super.setUp();
    mNano = new NanoRegression(SimilaritySvdTest.class);
  }

  @Override
  public void tearDown() throws Exception {
    super.tearDown();
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  public void testInit() {
    SimilaritySvd svd = new SimilaritySvd(5);
    assertEquals(-1, svd.getSvdDimension());
    assertEquals(-1, svd.getSvdRowsLength());
    try {
      svd.getSvdValue(0, 0);
      fail("returned value when not set up");
    } catch (IllegalStateException e) {
      // expected
    }
    try {
      svd.getSvdName(0);
      fail("returned value when not set up");
    } catch (IllegalStateException e) {
      // expected
    }
  }

  private SimilaritySvd readFromResource() throws IOException {
    SimilaritySvd svd = null;
    final String matrixString = FileHelper.resourceToString("com/rtg/similarity/resources/similarity.tsv");
    try (BufferedReader br = new BufferedReader(new StringReader(matrixString))) {
      String line;
      int count = 0;
      while ((line = br.readLine()) != null) {
        if (!line.startsWith("#")) {
          final String[] parts = line.split("\t");
          if (line.length() != 0) {
            if (count == 0) {
              //System.err.println(parts.length + " : " + parts[0]);
              svd = new SimilaritySvd(parts.length - 1);
            } else {
              svd.putName(count - 1, parts[0]);
              for (int i = 1; i < parts.length; i++) {
                svd.put(i - 1, count - 1, Integer.parseInt(parts[i]));
              }
            }
            count++;
          }
        }
      }
    }
    return svd;
  }

  public void testSerious() throws Exception {
    final File dir = FileUtils.createTempDir("similaritysvd", "test");
    try {
      SimilaritySvd svd = readFromResource();

      svd.decompose(3);
      assertEquals(3, svd.getSvdDimension());
      assertEquals(30, svd.getSvdRowsLength());

      final StringBuilder sb = new StringBuilder();
      for (int j = 0; j < svd.getSvdRowsLength(); j++) {
        for (int i = 0; i < 3; i++) {
          sb.append(svd.getSvdValue(j, i)).append(StringUtils.TAB);
        }
        sb.append(svd.getSvdName(j)).append(StringUtils.LS);
      }
      mNano.check("svd.tsv", sb.toString(), false);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testZeros() throws Exception {
    final long seed = new GregorianCalendar().getTimeInMillis() % 234567891L;
    final Random random = new Random(seed);
    final SimilaritySvd svd = readFromResource();

    int expected = 30;
    int[] rowsleft = new int[30];
    for (int i = 0; i < rowsleft.length; i++) {
      rowsleft[i] = i;
    }
    for (int i = 0; i < 20; i++) {
      final int r = random.nextInt(expected);
      final int row = rowsleft[r];
      rowsleft[r] = rowsleft[expected - 1];

      for (int j = 0; j < 30; j++) {
        svd.put(row, j, 0);
        svd.put(j, row, 0);
      }
      expected--;

      svd.decompose(3);
      assertEquals("seed=" + seed + " row=" + row, 3, svd.getSvdDimension());
      assertEquals("seed=" + seed + " row=" + row, expected, svd.getSvdRowsLength());

    }
  }

  public void testConstant() throws Exception {
    final long seed = new GregorianCalendar().getTimeInMillis() % 123456789L;
    final Random random = new Random(seed);
    final SimilaritySvd svd = readFromResource();

    int expected = 30;
    int[] rowsleft = new int[30];
    for (int i = 0; i < rowsleft.length; i++) {
      rowsleft[i] = i;
    }
    for (int i = 0; i < 20; i++) {
      final int r = random.nextInt(expected);
      final int row = rowsleft[r];
      rowsleft[r] = rowsleft[expected - 1];

      for (int j = 0; j < 30; j++) {
        svd.put(row, j, 42);
        svd.put(j, row, 42);
      }
      expected--;

      svd.decompose(3);
      assertEquals("seed=" + seed + " row=" + row, 3, svd.getSvdDimension());
      assertEquals("seed=" + seed + " row=" + row, expected, svd.getSvdRowsLength());
    }
  }
}
