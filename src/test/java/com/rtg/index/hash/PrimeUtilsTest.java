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

