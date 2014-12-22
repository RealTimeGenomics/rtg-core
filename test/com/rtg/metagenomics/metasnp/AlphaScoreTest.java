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
package com.rtg.metagenomics.metasnp;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 */
public class AlphaScoreTest extends TestCase {
  public void test() {
    final AlphaScore score = new AlphaScore(0.5, 0.2, 1, 2, 2, 2);
    assertTrue(Arrays.equals(new int[] {1, 2, 2, 2}, score.mCalls));
    assertEquals(0.5, score.mLikelihood, 0.0001);
    assertEquals(0.2, score.mBestPoss, 0.0001);
  }
}
