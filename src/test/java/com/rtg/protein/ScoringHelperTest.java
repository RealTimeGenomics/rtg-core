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
package com.rtg.protein;

import java.io.IOException;

import com.rtg.alignment.ActionsHelper;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.util.InvalidParamsException;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ScoringHelperTest extends TestCase {

  /**
   */
  public ScoringHelperTest(String name) {
    super(name);
  }

  public void testEScore() throws IOException, InvalidParamsException {
    final ProteinScoringMatrix matrix = new ProteinScoringMatrix();
    assertEquals(1.0402e-8, ScoringHelper.computeEScore(-100, 11, 100 * 1000, matrix), 0.0001e-8);
    assertEquals(0.041864, ScoringHelper.computeEScore(-67, 11, 60 * 1000 * 1000, matrix), 0.000001);
    assertEquals(5475.8724, ScoringHelper.computeEScore(-10, 21, 1000 * 1000, matrix), 0.001);
    assertEquals(115551.295, ScoringHelper.computeEScore(-3, 12, 1000 * 1000, matrix), 0.01);
    assertEquals(451000.000, ScoringHelper.computeEScore(0, 11, 1000 * 1000, matrix), 0.01);
  }

  public void testBitScore() throws IOException, InvalidParamsException {
    final ProteinScoringMatrix matrix = new ProteinScoringMatrix();
    assertEquals(30.4166, ScoringHelper.computeBitScore(-67, matrix), 0.0001);
    assertEquals(28.4906, ScoringHelper.computeBitScore(-62, matrix), 0.0001);
  }

  public void testIdentity() {
    final int[] sb = new int[20];
    assertEquals(0, ScoringHelper.percentId(sb));
    sb[ActionsHelper.ACTIONS_START_INDEX] = ActionsHelper.SAME << (ActionsHelper.ACTIONS_PER_INT - 1) * ActionsHelper.BITS_PER_ACTION;
    sb[ActionsHelper.ACTIONS_LENGTH_INDEX] = 1;
    assertEquals(100, ScoringHelper.percentId(sb));
    sb[ActionsHelper.ACTIONS_START_INDEX] |= ActionsHelper.MISMATCH << (ActionsHelper.ACTIONS_PER_INT - 2) * ActionsHelper.BITS_PER_ACTION;
    sb[ActionsHelper.ACTIONS_LENGTH_INDEX] = 2;
    assertEquals(50, ScoringHelper.percentId(sb));
    //next line is an identity operation that findbugs complains about
    //it is techincally correct to do the operation, but findbugs complains.
    //sb[ActionsHelper.ACTIONS_START_INDEX] |= ActionsHelper.SAME << 26;
    sb[ActionsHelper.ACTIONS_LENGTH_INDEX] = 3;
    assertEquals(67, ScoringHelper.percentId(sb));
    sb[ActionsHelper.ACTIONS_START_INDEX] |= ActionsHelper.MISMATCH << (ActionsHelper.ACTIONS_PER_INT - 4) * ActionsHelper.BITS_PER_ACTION;
    sb[ActionsHelper.ACTIONS_START_INDEX] |= ActionsHelper.MISMATCH << (ActionsHelper.ACTIONS_PER_INT - 5) * ActionsHelper.BITS_PER_ACTION;
    sb[ActionsHelper.ACTIONS_LENGTH_INDEX] = 5;
    assertEquals(40, ScoringHelper.percentId(sb));
    sb[ActionsHelper.ACTIONS_START_INDEX] |= ActionsHelper.MISMATCH << (ActionsHelper.ACTIONS_PER_INT - 6) * ActionsHelper.BITS_PER_ACTION;
    sb[ActionsHelper.ACTIONS_LENGTH_INDEX] = 6;
    assertEquals(33, ScoringHelper.percentId(sb));
  }

}
