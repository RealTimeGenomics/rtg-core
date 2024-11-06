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
