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

import com.rtg.mode.DnaUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

/**
 */
public class RcEditDistanceTest extends AbstractEditDistanceTest {

  @Override
  protected BidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty) {
    return new RcEditDistance(new GotohEditDistance(gapOpenPenalty, gapExtendPenalty, 1, 0, false));
  }

  public void testSomeLogging() {
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      Diagnostic.setLogStream(mps.printStream());
      try {
        final BidirectionalEditDistance ed = getEditDistanceInstance(1, 1);

        final byte[] read = new byte[]{1};
        final byte[] tmpl = new byte[]{2};


        ed.calculateEditDistance(read, read.length, tmpl, -1, false, 4, 4, true);
        ed.calculateEditDistance(read, read.length, tmpl, 2, false, 3, 4, true);
        ed.calculateEditDistance(read, read.length, tmpl, 3, false, 3, 4, true);
        ed.calculateEditDistance(read, read.length, tmpl, 1, false, 3, 4, true);
        ed.calculateEditDistance(read, -1, tmpl, 1, false, 4, 4, true);

        ed.logStats();

        TestUtils.containsAll(mps.toString(),
                                     "PROBLEM: RcEditDistance problem: 1 2/1 3 A",
                                     "PROBLEM: RcEditDistance problem: 1 3/1 3 A",
                                     "PROBLEM: RcEditDistance problem: 1 1/1 3 A",
                                     "PROBLEM: RcEditDistance problem: -1 1/1 4 A",
                                     "Can not create DPM for parameters rlen=0 bLength=1 dpmLength=0");
        assertEquals(mps.toString(), 10, mps.toString().trim().split("\n").length);
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }

  public void testSecondED() {
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      Diagnostic.setLogStream(mps.printStream());
      try {
        final GotohEditDistance ged = new GotohEditDistance(1, 1, 1, 0, false);
        final BidirectionalEditDistance ed = new RcEditDistance(new LoggingOnlyEditDistance(), ged);

        final int[] actions = ed.calculateEditDistance(DnaUtils.encodeString("aa"), 2, DnaUtils.encodeString("ttt"), 1, true, 5, 5, true);
        assertEquals(1, actions[ActionsHelper.TEMPLATE_START_INDEX]);
        assertEquals(0, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
        ged.logStats();
        TestUtils.containsAll(mps.toString(), "GotohEditDistance maxbytes=", "0 = 1");
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }
}
