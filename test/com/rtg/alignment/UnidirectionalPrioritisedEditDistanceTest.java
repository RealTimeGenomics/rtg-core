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
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParamsBuilder;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.reader.PrereadType;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.MaxShiftUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;


/**
 */
public class UnidirectionalPrioritisedEditDistanceTest extends AbstractUnidirectionalEditDistanceTest {

  @Override
  protected UnidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty, int subsPenalty, int unknownsPenalty) {
    return  new UnidirectionalPrioritisedEditDistance(new GotohEditDistance(gapOpenPenalty, gapExtendPenalty, subsPenalty, unknownsPenalty, false));
  }

  private class MinIntEd implements UnidirectionalEditDistance {
    private final int[] mActions;

    /**
     * Creates a Prioritised edit distance which tries several options in order
     */
    MinIntEd() {
      mActions = new int[12];
      mActions[ActionsHelper.ALIGNMENT_SCORE_INDEX] = Integer.MIN_VALUE;
    }

    @Override
    public void logStats() { }

    @Override
    public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
      if (rlen == 0) {
        return mActions;
      }
      if (maxScore < 0) {
        fail();
      }
      return null;
    }

    @Override
    public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos,
        int templateEndPos, int maxScore, int maxShift) {
      if (read.length == 0) {
        return mActions;
      }
      return null;
    }

    @Override
    public int[] calculateEditDistanceFixedEnd(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateExpectedStartPos,
        int templateEndPos, int maxScore, int maxShift) {
      if (read.length == 0) {
        return mActions;
      }
      return null;
    }

    @Override
    public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int maxScore, int maxShift) {
      if (read.length == 0) {
        return mActions;
      }
      return null;
    }
  }

  public void testPriority() {
    UnidirectionalEditDistance ed = new UnidirectionalPrioritisedEditDistance(new MinIntEd());

    final byte[] s1 = "gggggattttt".getBytes();
    final byte[] s2 = "gggggttttt".getBytes();
    DnaUtils.encodeArray(s1);
    DnaUtils.encodeArray(s2);

    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      int[] actions = ed.calculateEditDistance(s1, s1.length, s2, 0, 5, MaxShiftUtils.calculateDefaultMaxShift(s1.length), true);
      assertEquals(Integer.MAX_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      assertEquals(12, actions.length);

      ed.logStats();

      TestUtils.containsAll(mps.toString(),
          "calls=1", "nulls=1");

      actions = ed.calculateEditDistanceFixedBoth(s1, 0, s1.length, s2, 0, s2.length, 5, MaxShiftUtils.calculateDefaultMaxShift(s1.length));
      assertEquals(Integer.MAX_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      actions = ed.calculateEditDistanceFixedEnd(s1, 0, s1.length, s2, 0, s2.length, 5, MaxShiftUtils.calculateDefaultMaxShift(s1.length));
      assertEquals(Integer.MAX_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      actions = ed.calculateEditDistanceFixedStart(s1, 0, s1.length, s2, 0, 5, MaxShiftUtils.calculateDefaultMaxShift(s1.length));
      assertEquals(Integer.MAX_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

      ed = new UnidirectionalPrioritisedEditDistance(new MinIntEd(), new GotohEditDistance(1, 1, 1, 1, false));
      actions = ed.calculateEditDistance(s1, s1.length, s2, 0, 2, MaxShiftUtils.calculateDefaultMaxShift(s1.length), true);
      assertEquals(2, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

      actions = ed.calculateEditDistance(s1, 0, s2, 0, 5, MaxShiftUtils.calculateDefaultMaxShift(s1.length), true);
      assertEquals(Integer.MIN_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

      actions = ed.calculateEditDistanceFixedBoth(s1, 0, s1.length, s2, 0, s2.length, 2, 5);
      assertEquals(2, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      actions = ed.calculateEditDistanceFixedBoth(new byte[] {}, 0, 0, s2, 0, s2.length, 5, 5);
      assertEquals(Integer.MIN_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

      actions = ed.calculateEditDistanceFixedEnd(s1, 0, s1.length, s2, 0, s2.length, 2, 7);
      assertEquals(2, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      actions = ed.calculateEditDistanceFixedEnd(new byte[] {}, 0, 0, s2, 0, s2.length, 5, 7);
      assertEquals(Integer.MIN_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

      actions = ed.calculateEditDistanceFixedStart(s1, 0, s1.length, s2, 0, 2, 10);
      assertEquals(2, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
      actions = ed.calculateEditDistanceFixedStart(new byte[] {}, 0, 0, s2, 0, 5, 10);
      assertEquals(Integer.MIN_VALUE, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

      ed.calculateEditDistance(s1, 0, s2, 0, Integer.MAX_VALUE, 5, true);

      ed.logStats();

      TestUtils.containsAll(mps.toString(), "avg.score=2");
    } finally {
      mps.close();
    }
  }

  public void testSomeLogging() {
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      Diagnostic.setLogStream(mps.printStream());
      final UnidirectionalEditDistance ed = getEditDistanceInstance(1, 1, 1, 0);

      final byte[] read = new byte[]{1};
      final byte[] tmpl = new byte[]{2};


      ed.calculateEditDistance(read, read.length, tmpl, -1, 4, 4, true);
      ed.calculateEditDistance(read, read.length, tmpl, 2, 3, 4, true);
      ed.calculateEditDistance(read, read.length, tmpl, 3, 3, 4, true);
      ed.calculateEditDistance(read, read.length, tmpl, 1, 3, 4, true);
      try {
        ed.calculateEditDistance(read, -1, tmpl, 1, 4, 4, true);
      } catch (final Exception e) {
        assertTrue(e instanceof ArrayIndexOutOfBoundsException);
      }

      ed.logStats();

      TestUtils.containsAll(mps.toString(),
                                   "PROBLEM: UnidirectionalED problem: 1 2/1 3 A",
                                   "PROBLEM: UnidirectionalED problem: 1 3/1 3 A",
                                   "PROBLEM: UnidirectionalED problem: 1 1/1 3 A",
                                   "PROBLEM: UnidirectionalED problem: -1 1/1 4 A",
                                   "calls=5");
      final String[] sarr = mps.toString().split(" ");
      for (final String s : sarr) {
        if (s.startsWith("s)=")) {
          final String[] bla = s.split("=");
          assertEquals(2, bla.length);
          final Double d = Double.parseDouble(bla[1]);
          assertNotNull(d);
          assertTrue(d > 0.0);
        }
      }
    }
  }

  //this test demonstrates hop step edit distance's poisonous influence on the alignment chain. By changing the max score threshold we are able to align with a score better by 3!
  public void testHopStepBadMaxScore() {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      final byte[] read = DnaUtils.encodeString("CGTTTGAACCGGGGAGGCTGAGGTTGCCGTGAGCCAAGATTGTGCCACTGCACTCTAGCCTGGGCTACAGGGCAAGACTCCATTAAAAAAAAAAAAAAAC");
      final byte[] tmpl = DnaUtils.encodeString("CGTTTGAACCTGGGAGGCGGAGGTTGCAGTGAGCCAAGATTGTGCCACTGCACTCTAGCCTGGGCTACAGGGCAAGACTCCATTAAAAAAAAAAAAAAAAAACC");

      final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1)
                                       .outputParams(new NgsOutputParamsBuilder().filterParams(new NgsFilterParams.NgsFilterParamsBuilder().matedMaxMismatches(new IntegerOrPercentage("10%")).create()).create())
                                       .create();
      final EditDistance ed = EditDistanceFactory.createEditDistance(params, PrereadType.UNKNOWN, 100, 100);

      //HopStep can find an alignment with score 5 (the large seed at the end of the read ties the position, forcing a del of length 3)
      int[] actions = ed.calculateEditDistance(read, read.length, tmpl, 0, false, 10, 8, false);
      assertNotNull(actions);
      assertEquals(7, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);

      //providing a lower maxScore stops HopStep's alignment passing the threshold, and it falls through to SeededAligner which aligns with score 4 (all mismatches)
      actions = ed.calculateEditDistance(read, read.length, tmpl, 0, false, 6, 8, false);
      assertNotNull(actions);
      assertEquals(4, actions[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    } finally {
      Diagnostic.setLogStream();
    }
  }
}
