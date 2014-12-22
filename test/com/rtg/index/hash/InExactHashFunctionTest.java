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
    for (int i = 0; i < codes.length - windowSize; i++) {
      for (int j = i; j < i + windowSize - 1; j++) {
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

