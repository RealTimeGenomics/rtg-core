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

import java.io.IOException;

import com.rtg.mode.DnaUtils;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.MachineErrorParamsBuilder;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;

/**
 */
public class RcEditDistanceTest extends AbstractEditDistanceTest {

  @Override
  protected EditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty) {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(gapOpenPenalty)
        .gapExtendPenalty(gapExtendPenalty)
        .substitutionPenalty(1)
        .unknownsPenalty(0).create();
    return new RcEditDistance(new GotohEditDistance(params));
  }

  public void testAbstract() {
    final CgGotohEditDistanceTest test = new CgGotohEditDistanceTest();
    test.testCGOverlap0Rev();
    test.testCgOverlapRealWorld7();
  }

  public void testSomeLogging() {
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      Diagnostic.setLogStream(mps.printStream());
      try {
        final EditDistance ed = getEditDistanceInstance(1, 1);

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
        final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(0).create();
        final GotohEditDistance ged = new GotohEditDistance(params);
        final EditDistance ed = new RcEditDistance(new LoggingOnlyEditDistance(), ged);

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

  private static class TestGotohEditDistance extends CgGotohEditDistance {
    public TestGotohEditDistance(int maxShift, RealignParams params) {
      super(maxShift, params, 0);
    }
    @Override
    public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
      return ActionsHelper.build("=====B====X===============NNNNNN==========", 0, 2);
    }
  }

  public void testCGGotohHandling() {
    final RealignParamsImplementation mParams;
    try {
      mParams = new RealignParamsImplementation(new MachineErrorParamsBuilder().errors("cg_test_errors").create());
    } catch (final InvalidParamsException | IOException e) {
      throw new RuntimeException(e);
    }

    final RcEditDistance rced = new RcEditDistance(new TestGotohEditDistance(7, mParams));
    final int[] actions = rced.calculateEditDistance(DnaUtils.encodeString("attcttactanccccaactgagccc     agtaggagta".replaceAll(" ", "")), 35, DnaUtils.encodeString("atactcctacttttgctgggctcagttgggggaagtagaat"), 0, true, 25, 7, true);

    assertEquals(1, ActionsHelper.zeroBasedTemplateStart(actions));

  }
}
