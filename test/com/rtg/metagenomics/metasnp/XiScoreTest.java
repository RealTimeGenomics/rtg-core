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
public class XiScoreTest extends TestCase {
  public void test() {
    final XiScore score = new XiScore(0.5, 1.0, 2.0, 2.0, 2.0);
    assertTrue(Arrays.equals(new double[] {1, 2, 2, 2}, score.mXi));
    assertEquals(0.5, score.mScore, 0.0001);
  }
}
