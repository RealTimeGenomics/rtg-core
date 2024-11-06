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

package com.rtg.variant;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 */
public class Interpolate2dArrayTest extends TestCase {

  private static final int[][] CURVES = new int[][]{{-1, 4, -1, 6, -1, 10, -1},
    {-1, -1, -1, -1, -1, -1, -1},
    {-1, 20, -1, 15, -1, 9, -1},
  };
  private static final int[][] EXPECTED = new int[][]{{-1, 4, 5, 6, 8, 10, -1},
    {-1, -1, -1, -1, -1, -1, -1},
    {-1, 20, 18, 15, 12, 9, -1}
  };

  int[][] transpose(int[][] arr) {
    final int[][] result = new int[arr[0].length][arr.length];
    for (int i = 0; i < arr.length; ++i) {
      for (int j = 0; j < arr[i].length; ++j) {
        result[j][i] = arr[i][j];
      }
    }
    return result;
  }

  public void testColumn() {
    // These are expressed as rows because that's far easier to type/read
    final int[][] curves = transpose(CURVES);

    final int[][] expected = transpose(EXPECTED);


    Interpolate2dArray.column(curves, 0).interpolate();
    Interpolate2dArray.column(curves, 1).interpolate();
    Interpolate2dArray.column(curves, 2).interpolate();
    final String exp = Arrays.deepToString(expected);
    final String actual = Arrays.deepToString(curves);
    assertTrue(String.format("Expected <%s> but was <%s>", exp, actual), Arrays.deepEquals(expected, curves));
  }

  public void testRow() {
    // These are expressed as rows because that's far easier to type/read

    Interpolate2dArray.row(CURVES, 0).interpolate();
    Interpolate2dArray.row(CURVES, 1).interpolate();
    Interpolate2dArray.row(CURVES, 2).interpolate();
    final String exp = Arrays.deepToString(EXPECTED);
    final String actual = Arrays.deepToString(CURVES);
    assertTrue(String.format("Expected <%s> but was <%s>", exp, actual), Arrays.deepEquals(EXPECTED, CURVES));
  }

  public void testFill() {


  }

}
