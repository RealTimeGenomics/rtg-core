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

package com.rtg.sam.probe;

import htsjdk.samtools.CigarOperator;
import junit.framework.TestCase;

public class PositionAndStrandCheckerTest extends TestCase {

  public void testStats() {
    final PosChecker pos = new PosChecker(10);
    assertEquals(10, pos.mTolerance);
    assertEquals(21, pos.mPosDiffStats.length);
    pos.updateStrippedStats(CigarOperator.D, 7);
    pos.updateStrippedStats(CigarOperator.D, 7);
    pos.updateStrippedStats(CigarOperator.I, 5);
    pos.updateStrippedStats(CigarOperator.S, 9);
    pos.updateStrippedStats(CigarOperator.X, 13);
    pos.updateStrippedStats(CigarOperator.D, 30);

    assertEquals(1, pos.mMismatchStats[12]);
    assertEquals(2, pos.mDeletionStats[6]);
    assertEquals(1, pos.mInsertStats[4]);
    assertEquals(1, pos.mSoftClipStats[8]);
    assertEquals(1, pos.mDeletionStats[19]);
  }
}