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

import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractHashFunctionTest extends TestCase {

  protected abstract HashFunction getHashFunction(final int windowSize, final int bits);

  /**
   * Test method for {@link com.rtg.index.hash.ExactHashFunction}.
   */
  public final void test() {
    check(1, 1, new byte[] {0, 1});
    check(1, 2, new byte[] {0, 1, 2, 3});
    check(1, 5, new byte[] {0, 1, 17, 23, 31});
    check(1, 64, new byte[] {0, 1, 2, 3, 127, 126});
    check(2, 1, new byte[] {0, 1, 1, 0, 0});
    check(2, 32, new byte[] {0, 1, 7, 15, 31, 63, 127});
    check(12, 5, new byte[] {0, 1, 21, 24, 31, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 });
    check(12, 5, new byte[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
    check(32, 1, new byte[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, });
    check(32, 2, new byte[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, });
    check(64, 1, new byte[] {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, });
  }

  /**
   * Generate all the codes and check that they are all different.
   */
  protected void check(final int windowSize, final int bits, final byte[] codes) {
    final HashFunction hf = getHashFunction(windowSize, bits);
    Exam.integrity(hf);
    final Set<Long> al = new HashSet<>();
    for (int i = 0; i < windowSize; i++) {
      hf.hashStep(codes[i]);
    }
    assertTrue(al.add((long) hf.hashCode()));
    for (int i = windowSize; i < codes.length; i++) {
      assertTrue(al.add(hf.hashStep(codes[i])));
    }
  }
}

