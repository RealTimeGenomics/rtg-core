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

package com.rtg.variant;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 * @author kurt
 */
public class Interpolate2dArrayColumnTest extends TestCase {
  int[][] transpose(int[][] arr) {
    final int[][] result = new int[arr[0].length][arr.length];
    for (int i = 0; i < arr.length; i++) {
      for (int j = 0; j < arr[i].length; j++) {
        result[j][i] = arr[i][j];
      }
    }
    return result;
  }

  public void testSimple() {
    // These are expressed as rows because that's far easier to type/read
    int[][] curves = {{-1, 4, -1, 6, -1, 10, -1},
      {-1, -1, -1, -1, -1, -1, -1},
      {-1, 20, -1, 15, -1, 9, -1},
    };
    curves = transpose(curves);

    int[][] expected = {{-1, 4, 5, 6, 8, 10, -1},
      {-1, -1, -1, -1, -1, -1, -1},
      {-1, 20, 18, 15, 12, 9, -1}
    };
    expected = transpose(expected);


    new Interpolate2dArrayColumn(0, curves).process();
    new Interpolate2dArrayColumn(1, curves).process();
    new Interpolate2dArrayColumn(2, curves).process();
    final String exp = Arrays.deepToString(expected);
    final String actual = Arrays.deepToString(curves);
    assertTrue(String.format("Expected <%s> but was <%s>", exp, actual), Arrays.deepEquals(expected, curves));
  }

  public void testFill() {

  }

}