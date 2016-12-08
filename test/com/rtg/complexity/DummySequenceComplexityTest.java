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
package com.rtg.complexity;

import junit.framework.TestCase;

/**
 * Tests for AbstractSequenceComplexity class.
 *
 */
public class DummySequenceComplexityTest extends TestCase {

  private AbstractSequenceComplexity create(final int length) {
    return new AbstractSequenceComplexity(length) {
      @Override
      protected byte[] encodeString(final String str) {
        final byte[] b = new byte[str.length()];
        for (int i = 0; i < str.length(); ++i) {
          b[i] = (byte) str.charAt(i);
        }
        return b;
      }

      @Override
      protected double complexity(final byte[] b, int offset) {
        final int[] c = new int[26];
        for (int i = offset; i < b.length && i - offset < regionLength(); ++i) {
          c[b[i] - (byte) 'a']++;
        }
        double r = 0.0;
        for (int aC : c) {
          if (aC != 0) {
            r += 1.0 / aC;
          }
        }
        return r;
      }
    };
  }

  public void testConstrutor() {
    try {
      create(0);
      fail("Expected IAE");
    } catch (IllegalArgumentException iae) {
      // expected
    }
    try {
      create(21);
      fail("Expected IAE");
    } catch (IllegalArgumentException iae) {
      // expected
    }
    for (int i = 1; i <= 20; ++i) {
      create(i);
    }
  }

  public void testComplexity() {
    final AbstractSequenceComplexity sc = create(12);
    assertEquals(1.0 / 12.0, sc.minComplexity("aaaaaaaaaaaaaaaaaaaaa"));
    assertEquals(4.0 / 3.0, sc.minComplexity("acgtacgtacgt"));
    assertEquals(1.0 / 12.0, sc.maxComplexity("aaaaaaaaaaaaaaaaaaaaa"));
    assertEquals(4.0 / 3.0, sc.maxComplexity("acgtacgtacgt"));
  }
}
