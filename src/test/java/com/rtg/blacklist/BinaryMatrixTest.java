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
