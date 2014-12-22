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
package com.rtg.protein;

import junit.framework.TestCase;


/**
 */
public class ProteinBitScorerTest extends TestCase {

  private byte[] mRead;
  private byte[] mTmpl;

  @Override
  public void setUp() {
    mRead = new byte[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    mTmpl = new byte[] {1, 2, 3, 4, 5, 6, 1, 8, 7, 8, 9, 10};
  }

  @Override
  public void tearDown() {
    mRead = null;
    mTmpl = null;
  }

  public void testShiftScores1() {
    final int maxIndelLen = 1;
    final ProteinBitScorer bs = new ProteinBitScorer(maxIndelLen);
    assertEquals(1, bs.getMaxIndelLength());
    final int[] scores = bs.calculateFastScore(mRead, 0, 6, mTmpl, 0);
    assertEquals(2 * maxIndelLen + 2, scores.length);
    assertEquals(0, scores[maxIndelLen - 1]);
    assertEquals(6, scores[maxIndelLen]);
    assertEquals(0, scores[maxIndelLen + 1]);
  }

  public void testShiftScores2() {
    final int maxIndelLen = 2;
    final ProteinBitScorer bs = new ProteinBitScorer(maxIndelLen);
    assertEquals(2, bs.getMaxIndelLength());
    final int[] scores = bs.calculateFastScore(mRead, 0, 8, mTmpl, 0);
    assertEquals(2 * maxIndelLen + 2, scores.length);
    assertEquals(2, scores[maxIndelLen - 2]);
    assertEquals(0, scores[maxIndelLen - 1]);
    assertEquals(7, scores[maxIndelLen]);
    assertEquals(0, scores[maxIndelLen + 1]);
    assertEquals(0, scores[maxIndelLen + 2]);
  }

  public void testFastScore1() {
    final ProteinBitScorer bs = new ProteinBitScorer(1);

    checkTotalScore(6, bs.calculateFastScore(mRead, 0, 6, mTmpl, 0));
    checkTotalScore(6, bs.calculateFastScore(mRead, 0, 7, mTmpl, 0));
    checkTotalScore(7, bs.calculateFastScore(mRead, 0, 8, mTmpl, 0));
    checkTotalScore(7, bs.calculateFastScore(mRead, 0, 10, mTmpl, 0));

    checkTotalScore(0, bs.calculateFastScore(mRead, 0, 6, mTmpl, 2));
    checkTotalScore(4, bs.calculateFastScore(mRead, 0, 10, mTmpl, 2));
  }

  public void testFastScore2() {
    final ProteinBitScorer bs = new ProteinBitScorer(2);
    checkTotalScore(6, bs.calculateFastScore(mRead, 0, 6, mTmpl, 0));
    checkTotalScore(7, bs.calculateFastScore(mRead, 0, 7, mTmpl, 0));
    checkTotalScore(8, bs.calculateFastScore(mRead, 0, 8, mTmpl, 0));
    checkTotalScore(10, bs.calculateFastScore(mRead, 0, 10, mTmpl, 0));
  }

  private void checkTotalScore(final int expected, final int[] scores) {
    assertEquals(expected, scores[scores.length - 1]);
  }

  public void testFillup() {
    final ProteinBitScorer bs = new ProteinBitScorer(0);
    final byte[] data = {2, 31, 15, 7, 3, 1};
    final long[] result = new long[5];
    bs.fillup(data, 1, 6, 0, result);
    assertEquals(0x3E, result[0]);
    assertEquals(0x3C, result[1]);
    assertEquals(0x38, result[2]);
    assertEquals(0x30, result[3]);
    assertEquals(0x20, result[4]);
  }
}
