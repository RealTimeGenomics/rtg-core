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

import java.util.HashSet;
import java.util.Set;


/**
 */
public class InExactHashFunctionTest extends AbstractHashFunctionTest {
  private static final byte[] LONG_CODE = //128 nt
    new byte[] {
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 3,
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
  };

  private static final byte[] LONG_CODE_REPEAT = //128 nt
    new byte[] {
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
    0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2, 0, 1, 3, 2,
    0, 1, 3,
  };

  @Override
  protected HashFunction getHashFunction(final int windowSize, final int bits) {
    return new InExactHashFunction(windowSize);
  }

  /**
   * Test method for {@link com.rtg.index.hash.ExactHashFunction}.
   */
  public final void testInexact() {
    check(64, 1, LONG_CODE);
    check(64, 2, LONG_CODE);
    check(64, 2, LONG_CODE_REPEAT);
  }

  @Override
  protected void check(final int windowSize, final int bits, final byte[] codes) {
    final InExactHashFunction hf = (InExactHashFunction) getHashFunction(windowSize, bits);
    hf.integrity();
    final Set<Long> al = new HashSet<>();
    for (int i = 0; i < codes.length - windowSize; ++i) {
      for (int j = i; j < i + windowSize - 1; ++j) {
        hf.hashStep(codes[j]);
        hf.integrity();
      }
      assertTrue("" + i, al.add(hf.hashStep(codes[i + windowSize - 1])));
      hf.integrity();
      hf.reset();
      hf.integrity();
    }
  }

  /**
   * Test method for {@link com.rtg.util.integrity.IntegralAbstract#toString()}.
   */
  public final void testToString() {
    final InExactHashFunction hf = (InExactHashFunction) getHashFunction(3, 2);
    hf.integrity();
    assertEquals("IrvineHashFunction window size=3 hash=0 i=0", hf.toString());
    hf.hashStep((byte) 3);
    final String expectedString = "IrvineHashFunction window size=3 hash=6137546356583794141 i=1";
    assertEquals(expectedString, hf.toString());
  }

}

