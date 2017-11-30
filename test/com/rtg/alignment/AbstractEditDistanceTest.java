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
package com.rtg.alignment;

import com.rtg.AbstractTest;
import com.rtg.mode.DnaUtils;

/**
 */
public abstract class AbstractEditDistanceTest extends AbstractTest {

  protected abstract BidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty);

  public void testStartPositionSet() throws Exception {
    final byte[] read = DnaUtils.encodeString("acgta");
    final byte[] template = DnaUtils.encodeString("tttttacgtatc");

    final BidirectionalEditDistance ed = getEditDistanceInstance(1, 1);

    int[] actions = ed.calculateEditDistance(read, read.length, template, 5, false, 1, 1, false);

    if (actions != null) {
      assertEquals(0, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      assertEquals(5, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }

    actions = ed.calculateEditDistance(read, read.length, template, 6, false, 1, 1, false);

    if (actions != null) {
      assertEquals(0, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      assertEquals(5, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }

    actions = ed.calculateEditDistance(read, read.length, template, 7, false, 1, 1, false);

    if (actions != null) {
      assertEquals(Integer.MAX_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      assertEquals(7, actions[ActionsHelper.TEMPLATE_START_INDEX]);
    }
  }
}
