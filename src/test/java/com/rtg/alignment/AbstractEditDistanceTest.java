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
package com.rtg.alignment;

import com.rtg.AbstractTest;
import com.rtg.mode.DnaUtils;

/**
 */
public abstract class AbstractEditDistanceTest extends AbstractTest {

  protected abstract BidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty);

  public void testStartPositionSet() {
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
