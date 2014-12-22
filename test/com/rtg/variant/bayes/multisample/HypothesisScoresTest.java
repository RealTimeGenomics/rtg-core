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

import junit.framework.TestCase;

/**
 */
public class HypothesisScoresTest extends TestCase {

  public void test() {
    final HypothesisScores hs = new HypothesisScores(null, false, 0.42, null);
    assertEquals(0.42, hs.getNonIdentityPosterior());
    assertFalse(hs.isInteresting());
    assertNull(hs.getBs());
  }
}
