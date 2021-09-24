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

package com.rtg.variant.bayes.multisample;

import com.rtg.variant.bayes.MockGenotypeMeasure;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

import junit.framework.TestCase;

/**
 */
public class HypothesisScoreTest extends TestCase {

  public void test() {
    final HypothesisScore hs = new HypothesisScore(new MockGenotypeMeasure(0, 42, LogApproximatePossibility.SINGLETON.ln2Poss(0.42), Double.NaN));
    assertEquals(42, hs.hypothesis());
    assertEquals(0.42, hs.posterior());
    assertTrue(Double.isNaN(hs.nonIdentityPosterior()));
  }
}
