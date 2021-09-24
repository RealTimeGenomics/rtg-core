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

package com.rtg.metagenomics;

import junit.framework.TestCase;

/**
 */
public class LineTest extends TestCase {

  private static class LineTester extends Line {
    @Override
    public double value(double delta) {
      return 3.0 * delta + 42.0;
    }
  }

  public void test() {
    final LineTester lt = new LineTester();
    assertEquals(0, lt.derivativeOrder());
    assertEquals(42.0, lt.value(0.0));
    final double[] ltd = lt.values(0.0);
    assertEquals(1, ltd.length);
    assertEquals(42.0, ltd[0]);
  }

}
