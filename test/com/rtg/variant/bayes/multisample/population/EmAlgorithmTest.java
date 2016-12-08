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

package com.rtg.variant.bayes.multisample.population;

import com.rtg.variant.bayes.MockGenotypeMeasure;
import com.rtg.variant.bayes.multisample.HypothesisScore;

import junit.framework.TestCase;

/**
 */
public class EmAlgorithmTest extends TestCase {

  public void testCallToStringSmall1() {
    final HypothesisScore[] calls = new HypothesisScore[15];
    for (int i = 0; i < calls.length; ++i) {
      calls[i] = new HypothesisScore(new MockGenotypeMeasure(0, i % 10, 0.0, 0));
    }
    final String str = EmAlgorithm.callToString(10, calls);
    assertEquals("0123456789 01234", str);
  }

  //long enough to get two complete cycles -- tests % operation
  public void testCallToStringSmall2() {
    final HypothesisScore[] calls = new HypothesisScore[21];
    for (int i = 0; i < calls.length; ++i) {
      calls[i] = new HypothesisScore(new MockGenotypeMeasure(0, i % 10, 0.0, 0));
    }
    final String str = EmAlgorithm.callToString(10, calls);
    assertEquals("0123456789 0123456789 0", str);
  }

  public void testCallToStringBig() {
    final HypothesisScore[] calls = new HypothesisScore[15];
    for (int i = 0; i < calls.length; ++i) {
      calls[i] = new HypothesisScore(new MockGenotypeMeasure(0, i, 0.0, 0));
    }
    final String str = EmAlgorithm.callToString(11, calls);
    assertEquals("  0  1  2  3  4  5  6  7  8  9 |  10 11 12 13 14", str);
  }
}
