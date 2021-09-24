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

package com.rtg.blacklist;

import java.util.Random;

import junit.framework.TestCase;

public class BinaryMatrixTest extends TestCase {


  public void test() {
    final Random r = new Random(42);
    for (int k = 1; k < 64; ++k) {
      final BinaryMatrix matrix = BinaryMatrix.createReversibleMatrix(k);
      final BinaryMatrix inverse = matrix.invert();

      final long maxVector = (1L << k) - 1;
      final long unVector = ~maxVector;

      for (int i = 0; i < 20; ++i) {
        final long val = random(r, k);
        final long shuffled = matrix.times(val);
        assertEquals(0, shuffled & unVector);
        assertEquals(val, inverse.times(shuffled));
      }
    }
  }

  static long random(Random r, int numBits) {
    long ret = 0;
    for (int i = 0; i < numBits; ++i) {
      ret |= (r.nextBoolean() ? 1L : 0L) << i;
    }
    return ret;
  }
}
