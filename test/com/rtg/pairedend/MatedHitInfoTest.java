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
package com.rtg.pairedend;

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.AlignmentResult;

import junit.framework.TestCase;

/**
 */
public class MatedHitInfoTest extends TestCase {

  public void test() {
    final MatedHitInfo info = new MatedHitInfo();
    assertNull(info.mAlignmentLeft);
    assertNull(info.mAlignmentRight);
    assertEquals(-1, info.mAlignmentScoreLeft);
    assertEquals(-1, info.mAlignmentScoreRight);
    info.setValues(1, false, true, 3, false, 2);
    assertFalse(info.mFirstRight);
    assertFalse(info.mReverseComplementRight);
    assertTrue(info.mReverseComplementLeft);
    info.setValues(1, false, false, 3, true, 2);
    assertFalse(info.mFirstRight);
    assertTrue(info.mReverseComplementRight);
    assertFalse(info.mReverseComplementLeft);
    info.setValues(1, true, false, 3, false, 2);
    assertTrue(info.mFirstRight);
    assertFalse(info.mReverseComplementRight);
    assertFalse(info.mReverseComplementLeft);
    assertEquals(1L, info.mReadId);
    assertEquals(2, info.mTemplateStartRight);
    assertEquals(3, info.mTemplateStartLeft);
    assertNull(info.mAlignmentLeft);
    assertNull(info.mAlignmentRight);
    final int[] actionsl = ActionsHelper.build("", 0, 12);
    final int[] actionsr = ActionsHelper.build("", 0, 18);
    final AlignmentResult left = new AlignmentResult(null, actionsl, null);
    final AlignmentResult right = new AlignmentResult(null, actionsr, null);
    info.setAlignments(left, right);
    assertTrue(left == info.mAlignmentLeft);
    assertTrue(right == info.mAlignmentRight);
    assertEquals(12, info.mAlignmentScoreLeft);
    assertEquals(18, info.mAlignmentScoreRight);
    info.setValues(1, true, false, 3, false, 2);
    assertNull(info.mAlignmentLeft);
    assertNull(info.mAlignmentRight);
    assertEquals(-1, info.mAlignmentScoreLeft);
    assertEquals(-1, info.mAlignmentScoreRight);
    assertFalse(info.mLeftIsVeryPoor);
  }

}
