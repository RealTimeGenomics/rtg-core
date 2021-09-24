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
package com.rtg.position.output;

import junit.framework.TestCase;

/**
 */
public class ScannerTest extends TestCase {

  private static final Scanner<String> SCANNER = new Scanner<String>() {

    private StringBuilder mStringBuilder;

    @Override
    protected void init() {
      mStringBuilder = new StringBuilder();
    }
    @Override
    protected String result() {
      return mStringBuilder.toString();
    }

    @Override
    protected double resultScore() {
      return Double.NEGATIVE_INFINITY;
    }


    @Override
    protected void scan(final String region, final boolean forward, final long lo, final long hi) {
      mStringBuilder.append(lo).append("...").append(hi).append(":");
    }
  };

  //bug seen in the wild
  public void testBug1() {
    check(10, 20,  7, 15, "10...15:18...20:");
    check(10, 20, -9, -4, "13...18:");
    check(10, 20, -8, -3, "14...19:");
    check(10, 20, -7, -2, "15...20:");
    check(10, 20, -6, -1, "10...10:16...20:");
    check(10, 20, -5,  0, "10...11:17...20:");
    check(10, 20, -4,  4, "10...15:18...20:");
    check(10, 20, -3,  5, "10...16:19...20:");
    check(10, 20, -2,  6, "10...17:20...20:");
    check(10, 20, -1,  7, "10...18:");
  }

  //bug seen in the wild
  public void testBug2() {
    check(10, 20,  7, 15, "10...15:18...20:");
    check(10, 20, 21, 29, "10...18:");
    check(10, 20, 22, 30, "11...19:");
    check(10, 20, 23, 31, "12...20:");
    check(10, 20, 24, 32, "10...10:13...20:");
    check(10, 20, 25, 33, "10...11:14...20:");
  }

//bug seen in the wild
  public void testBug3() {
    check(5, 9,  -11, -11, "9...9:");
    check(5, 9,  -9, -9, "6...6:");
    check(5, 9,  -10, -10, "5...5:");
  }

  public void testBug4() {
    check(5, 9,  9, 9, "9...9:");
    check(5, 9,  10, 10, "5...5:");
    check(5, 9,  11, 11, "6...6:");
    check(5, 9,  12, 12, "7...7:");
    check(5, 9,  13, 13, "8...8:");
    check(5, 9,  14, 14, "9...9:");
    check(5, 9,  15, 15, "5...5:");
  }


  public void test1() {
    checkRange(7, 15, "10...15:18...20:");
    checkRange(8, 16, "10...16:19...20:");
    checkRange(9, 17, "10...17:20...20:");
    checkRange(10, 18, "10...18:");
    checkRange(11, 19, "11...19:");
    checkRange(12, 20, "12...20:");
    checkRange(13, 21, "10...10:13...20:");
    checkRange(14, 22, "10...11:14...20:");
    checkRange(15, 23, "10...12:15...20:");
    checkRange(16, 24, "10...13:16...20:");
    checkRange(17, 25, "10...14:17...20:");
    checkRange(18, 26, "10...15:18...20:");
  }

  private void checkRange(final int lo, final int hi, final String expected) {
    for (int i = -5; i <= 5; ++i) {
      final int j = 11 * i;
      check(10, 20, lo + j, hi + j, expected);
    }
  }

  public void test2() {
    check(10, 20, 10, 20, "10...20:");
    check(10, 20, 10, 20, "10...20:");
    check(10, 20, 10, 20, "10...20:");
    check(10, 20, 10, 19, "10...19:");
    check(10, 20, 11, 20, "11...20:");
    check(10, 20,  9, 20, "10...20:");
    check(10, 20, 12, 21, "10...10:12...20:");
    check(10, 20, 12, 22, "10...20:");
    check(10, 20, 12, 23, "10...20:");
    check(10, 20, 15, 20, "15...20:");
    check(10, 20, 15, 21, "10...10:15...20:");
    check(10, 20, 15, 22, "10...11:15...20:");
    check(10, 20, 15, 24, "10...13:15...20:");
    check(10, 20, 15, 25, "10...20:");
    check(10, 20, 15, 25, "10...20:");
    check(10, 20, 10, 15, "10...15:");
    check(10, 20,  9, 15, "10...15:20...20:");
    check(10, 20,  8, 15, "10...15:19...20:");
    check(10, 20,  6, 15, "10...15:17...20:");
    check(10, 20,  5, 15, "10...20:");
  }

  /* Check negative numbers */
  public void test3() {
    check(0, 10, 0, 10, "0...10:");
    check(0, 10, -1, 9, "0...10:");
    check(0, 10, -1, 5, "0...5:10...10:");
  }

  /* Check outside range */
  public void test4() {
    check(0, 13, -2, -1, "12...13:");
    check(0, 4, -1, 0, "0...0:4...4:");
  }

  private void check(final long first, final long last, final long lo, final long hi, final String exp) {
    assertEquals(exp, SCANNER.scanAll("region", true, first, last, lo, hi));
  }
}
