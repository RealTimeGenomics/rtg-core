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
