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
package com.rtg.index.hash;


/**
 */
//TODO lots more testing
public class RandomizedExactHashFunctionTest extends AbstractHashFunctionTest {


  @Override
  protected HashFunction getHashFunction(final int windowSize, final int bits) {
    return new RandomizedExactHashFunction(windowSize, bits);
  }

  /**
   * Test method for {@link com.rtg.index.hash.ExactHashFunction#reset()}.
   */
  public final void testReset() {
    final RandomizedExactHashFunction hf = (RandomizedExactHashFunction) getHashFunction(3, 2);
    hf.integrity();
    //see testReset in ExactHashFunctionTest
    assertEquals(0, hf.hashStep((byte) 0));
    assertTrue(1 != hf.hashStep((byte) 1));
    assertTrue(5 != hf.hashStep((byte) 1));
    assertTrue(21 != hf.hashStep((byte) 1));
    assertTrue(21 != hf.hashStep((byte) 1));
    assertTrue(22 != hf.hashStep((byte) 2));
    assertTrue(27 != hf.hashStep((byte) 3));
    hf.reset();
    hf.integrity();
    assertEquals(0, hf.hashStep((byte) 0));
    assertTrue(1 != hf.hashStep((byte) 1));
    assertEquals("RandomizedExactHashFunction window size=3 bits=2 mask=63 hash=1 prime=3", hf.toString());
  }

  /**
   * Check that the result of the hash correctly inverts for different numbers and each
   * possible bit length.
   */
  public void testHash() {
    for (int i = 1; i < 64; ++i) {
      check(i, 1L);
      check(i, (1L << i) - 1L);
      check(i, ((1L << i) - 1L) / 3);
    }
    check(64, 1L);
    check(64, -1L);
    check(64, Long.MAX_VALUE / 3);
  }

  /**
   * Check that constructing the specified hash one bit at a time and then
   * inverting the prime multiplier ends up with the original value.
   */
  private void check(final int bits, final long hash) {
    //System.err.println("check bits=" + bits + " hash=" + hash);
    final RandomizedExactHashFunction rh = new RandomizedExactHashFunction(bits, 1);
    long x = 0;
    for (int i = 0; i < bits; ++i) {
      x = rh.hashStep((byte) (hash >> (bits - i - 1) & 1L));
    }
    if (bits == 64) {
      final long y = x * PrimeUtils.primeInverse(bits);
      assertEquals("bits=" + bits + " hash=" + hash + " x=" + x, hash, y);
    } else {
      final long y = (x * PrimeUtils.primeInverse(bits)) & ((1L << bits) - 1L);
      assertEquals("bits=" + bits + " hash=" + hash + " x=" + x, hash, y);
    }
  }
}

