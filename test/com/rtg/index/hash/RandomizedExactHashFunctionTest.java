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
    for (int i = 1; i < 64; i++) {
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
    for (int i = 0; i < bits; i++) {
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

