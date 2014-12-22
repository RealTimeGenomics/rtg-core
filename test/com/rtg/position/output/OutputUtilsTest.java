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
package com.rtg.position.output;

import com.rtg.position.output.OutputUtils.ScoreResult;

import junit.framework.TestCase;

/**
 */
public class OutputUtilsTest extends TestCase {

  public final void testGapScore() {
    check(4, 4, 4, 4, 12.5, 3);
    check(8, 8, 4, 4, 25.0, 6);
    check(4, 3, 4, 4, 2.5, 3);
    check(3, 4, 4, 4, 2.5, 3);
    check(8, 7, 4, 4, 15.0, 6);
    check(7, 8, 4, 4, 15.0, 6);
    check(8, 9, 4, 4, 19.5, 7);
    check(9, 8, 4, 4, 19.5, 7);

    check(12, 12, 4, 4, 37.5, 9);
    check(12, 13, 4, 4, 32.0, 10);
    check(13, 12, 4, 4, 32.0, 10);

  }

  void check(final int aGap, final int bGap, final int wordSize, final int stepSize, final double expScore, final int expId) {
    final ScoreResult res = new ScoreResult();
    OutputUtils.gapScore(res, aGap, bGap, wordSize, stepSize);
    assertEquals(expScore, res.score());
    assertEquals(expId, res.idCount());
  }

}
