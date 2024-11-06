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

package com.rtg.similarity;


import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.GregorianCalendar;
import java.util.Random;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.util.StringUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 * Tests for SimilaritySvd class.
 */
public class SimilaritySvdTest extends AbstractNanoTest {

  public void testInit() {
    final SimilaritySvd svd = new SimilaritySvd(5);
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
              for (int i = 1; i < parts.length; ++i) {
                svd.put(i - 1, count - 1, Integer.parseInt(parts[i]));
              }
            }
            ++count;
          }
        }
      }
    }
    return svd;
  }

  public void testSerious() throws Exception {
    final File dir = FileUtils.createTempDir("similaritysvd", "test");
    try {
      final SimilaritySvd svd = readFromResource();

      svd.decompose(3);
      assertEquals(3, svd.getSvdDimension());
      assertEquals(30, svd.getSvdRowsLength());

      final StringBuilder sb = new StringBuilder();
      for (int j = 0; j < svd.getSvdRowsLength(); ++j) {
        for (int i = 0; i < 3; ++i) {
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
    final int[] rowsleft = new int[30];
    for (int i = 0; i < rowsleft.length; ++i) {
      rowsleft[i] = i;
    }
    for (int i = 0; i < 20; ++i) {
      final int r = random.nextInt(expected);
      final int row = rowsleft[r];
      rowsleft[r] = rowsleft[expected - 1];

      for (int j = 0; j < 30; ++j) {
        svd.put(row, j, 0);
        svd.put(j, row, 0);
      }
      --expected;

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
    final int[] rowsleft = new int[30];
    for (int i = 0; i < rowsleft.length; ++i) {
      rowsleft[i] = i;
    }
    for (int i = 0; i < 20; ++i) {
      final int r = random.nextInt(expected);
      final int row = rowsleft[r];
      rowsleft[r] = rowsleft[expected - 1];

      for (int j = 0; j < 30; ++j) {
        svd.put(row, j, 42);
        svd.put(j, row, 42);
      }
      --expected;

      svd.decompose(3);
      assertEquals("seed=" + seed + " row=" + row, 3, svd.getSvdDimension());
      assertEquals("seed=" + seed + " row=" + row, expected, svd.getSvdRowsLength());
    }
  }
}
