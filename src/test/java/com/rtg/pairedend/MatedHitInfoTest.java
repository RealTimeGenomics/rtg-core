/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
