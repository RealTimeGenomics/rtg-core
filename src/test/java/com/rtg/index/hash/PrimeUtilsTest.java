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

import junit.framework.TestCase;

/**
 */
public class PrimeUtilsTest extends TestCase {


  /** Make sure that each of the original primes is invertible. */
  public void testValid() {
    PrimeUtils.integrity();
    for (int i = 1; i < 64; ++i) {
      final long x = (PrimeUtils.prime(i) * PrimeUtils.primeInverse(i)) & (((1L << i) - 1L));
      assertEquals(i + ":" + x, 1L, x);
    }
    final long x2 = PrimeUtils.prime(64) * PrimeUtils.primeInverse(64);
    assertEquals(64 + ":" + x2, 1L, x2);
 }

  public void testBad() {
    try {
      PrimeUtils.prime(0);
      fail();
    } catch (final IllegalArgumentException e) {
        assertEquals("Number of bits should be between 1 and 64 inclusive. bits=0", e.getMessage());
    }

    try {
      PrimeUtils.primeInverse(0);
      fail();
    } catch (final IllegalArgumentException e) {
        assertEquals("Number of bits should be between 1 and 64 inclusive. bits=0", e.getMessage());
    }

    try {
      PrimeUtils.prime(65);
      fail();
    } catch (final IllegalArgumentException e) {
        assertEquals("Number of bits should be between 1 and 64 inclusive. bits=65", e.getMessage());
    }

    try {
      PrimeUtils.primeInverse(65);
      fail();
    } catch (final IllegalArgumentException e) {
        assertEquals("Number of bits should be between 1 and 64 inclusive. bits=65", e.getMessage());
    }
}

}

