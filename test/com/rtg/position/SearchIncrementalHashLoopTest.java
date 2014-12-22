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
package com.rtg.position;

import java.io.IOException;

import com.rtg.index.Index;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.index.hash.HashLoop;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderLongMock;

/**
 */
public class SearchIncrementalHashLoopTest extends AbstractPositionHashLoopTest {

  private static FinderPositionOutput getOutputVarsEmpty() {
    return new FinderPositionOutput(null, new MockPositionOutput());
  }
  @Override
  protected HashLoop getHashLoop1(final int windowSize, int stepSize,
      final int maxId, final int bits, final int[] count) {
    final HashLoop hashLoop;
    final int[] mi = new int[1];
    mi[0] = maxId;
    hashLoop =
      new SearchIncrementalHashLoop(stepSize, new ExactHashFunction(windowSize, bits), getOutputVarsEmpty(), (FinderPositionOutput) null, (Index) null, false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        assertTrue("" + internalId, internalId >= 0 && internalId < mi[0]);
        count[0]++;
      }

      @Override
      public void hashCallBidirectional(final long hashFoward, final long hashReverse, final int stepPosition, final int internalId) {
        hashCall(hashFoward, internalId, stepPosition);
        hashCall(hashReverse, internalId + 1, stepPosition);
      }

    };
    return hashLoop;
  }

  @Override
  protected HashLoop getHashLoop1a(final int windowSize, int stepSize, final int bits) {
    final HashLoop hashLoop1;
    hashLoop1 =
      new SearchIncrementalHashLoop(stepSize, new ExactHashFunction(windowSize, bits), getOutputVarsEmpty(), (FinderPositionOutput) null, (Index) null, false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        // do nothing
      }

      @Override
      public void hashCallBidirectional(final long hashFoward, final long hashReverse, final int stepPosition, final int internalId) {
      }


    };
    return hashLoop1;
  }

  @Override
  protected HashLoop getHashLoop3(final int windowSize, int stepSize,
      final long[] expectedL, final int[] expectedI, final int bits, final int[] count) {
    final HashLoop hashLoop;
    hashLoop =
      new SearchIncrementalHashLoop(stepSize, new ExactHashFunction(windowSize, bits, true), getOutputVarsEmpty(), (FinderPositionOutput) null, (Index) null, false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        assertEquals(count[0] + "", expectedL[count[0]], hash);
        assertEquals(count[0] + "", expectedI[count[0]], internalId);
        count[0]++;
      }

      @Override
      public void hashCallBidirectional(final long hashFoward, final long hashReverse, final int stepPosition, final int internalId) {
        hashCall(hashFoward, internalId, stepPosition);
        hashCall(hashReverse, internalId + 1, stepPosition);
      }


    };
    return hashLoop;
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   */
  @Override
  public final void testReadLength() {

  }

  private static final class TestLongIncrementalHashLoop extends SearchIncrementalHashLoop {
    public TestLongIncrementalHashLoop(int stepSize, final ExactHashFunction function) {
      super(stepSize, function, getOutputVarsEmpty(), null, null, false);
    }

    @Override
    public void hashCall(final long hash, final int internalId, final int stepPosition) {
    }
  }

  @Override
  protected void getLongLoop() throws IOException {
    final HashLoop hashLoop =
      new TestLongIncrementalHashLoop(1, new ExactHashFunction(1, 2));
    final ReaderParams re = new MockReaderParams(new ReaderLongMock(Integer.MAX_VALUE), SequenceMode.UNIDIRECTIONAL);
    final ISequenceParams se = new MockSequenceParams(re , 0, 0);
    hashLoop.execLoop(se, HashLoop.makeBuffer(se.reader()));
  }

  @Override
  public void testExceptionMessage() {

  }


  /*
  @Override
  public final void test4() {
    final String str = ""
      + ">x1" + StringUtils.LS
      + "ACGTACG" + StringUtils.LS;
    //0, 1, 2, 3, 0, 1, 2
    //3, 2, 1, 0, 3, 2, 1
    final SequencesReader sr = getReaderDNA(str);
    try {
      checkCount3(sr, 1, 1, UNIDIRECTIONAL,
          new long[] {0, 1, 2, 3, 0, 1, 2},
          new int[]  {0, 0, 0, 0, 0, 0, 0}
          );
      checkCount3(sr, 1, 1, BIDIRECTIONAL,
          //new long[] {0, 1, 2, 3, 0, 1, 2,
          //            1, 2, 3, 0, 1, 2, 3},
          new long[] {0, 3, 1, 2, 2, 1, 3, 0, 0, 3, 1, 2, 2, 1},
          new int[]  {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1}
          );
      checkCount3(sr, 1, 1, TRANSLATED,
          new long[] {16, 18, 1, 16, 19, 1, 16, 19, 1, 18},
          new int[]  {0,  0, 1,  1,  2, 3,  3,  4, 4,  5}
          );
      checkCount3(sr, 2, 1, UNIDIRECTIONAL,
          new long[] {1, 6, 11, 12, 1, 6},
          new int[]  {0, 0,  0,  0, 0, 0}
          );
      checkCount3(sr, 2, 1, BIDIRECTIONAL,
          //new long[] {1, 6, 11, 12, 1, 6,
          //            6, 11, 12, 1, 6, 11},
          //          AC|GT|CG|CG| GT|AC| TA| TA|AC|GT| CG| CG
          new long[] {1, 11, 6, 6, 11, 1, 12, 12, 1, 11, 6,  6},
          new int[]  {0, 1,  0,  1, 0, 1, 0,  1,  0, 1, 0,  1}
          );
      checkCount3(sr, 2, 1, TRANSLATED,
          new long[] {530, 48, 48, 609},
          new int[]  {0,   1,  3,   4}
          );
    } finally {
      sr.close();
    }
  }

  public final void test5() {
    Diagnostic.setLogStream();
    final String str = ""
      + ">x" + StringUtils.LS
      + "ACGT" + StringUtils.LS;
    final SequencesReader sr = getReaderDNA(str);
    try {
      checkCount3(sr, 2, 2, SequenceMode.BIDIRECTIONAL,
          new long[] {1, 11, 11, 1},
          new int[]  {0, 1,  0, 1}
          );
    } finally {
      sr.close();
    }
  }
   */
}
