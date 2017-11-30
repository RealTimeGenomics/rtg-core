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
