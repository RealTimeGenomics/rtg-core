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
  protected HashLoop getHashLoop1(final int windowSize, int stepSize, final int maxId, final int bits, final int[] count) {
    final int[] mi = new int[1];
    mi[0] = maxId;
    return new SearchIncrementalHashLoop(stepSize, new ExactHashFunction(windowSize, bits), getOutputVarsEmpty(), null, null, false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        assertTrue("" + internalId, internalId >= 0 && internalId < mi[0]);
        ++count[0];
      }

      @Override
      public void hashCallBidirectional(final long hashForward, final long hashReverse, final int stepPosition, final int internalId) {
        hashCall(hashForward, internalId, stepPosition);
        hashCall(hashReverse, internalId + 1, stepPosition);
      }
    };
  }

  @Override
  protected HashLoop getHashLoop1a(final int windowSize, int stepSize, final int bits) {
    return new SearchIncrementalHashLoop(stepSize, new ExactHashFunction(windowSize, bits), getOutputVarsEmpty(), null, null, false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        // do nothing
      }

      @Override
      public void hashCallBidirectional(final long hashForward, final long hashReverse, final int stepPosition, final int internalId) {
      }
    };
  }

  @Override
  protected HashLoop getHashLoop3(final int windowSize, int stepSize,
                                  final long[] expectedL, final int[] expectedI, final int bits, final int[] count) {
    return new SearchIncrementalHashLoop(stepSize, new ExactHashFunction(windowSize, bits, true), getOutputVarsEmpty(), null, null, false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        assertEquals(count[0] + "", expectedL[count[0]], hash);
        assertEquals(count[0] + "", expectedI[count[0]], internalId);
        ++count[0];
      }

      @Override
      public void hashCallBidirectional(final long hashForward, final long hashReverse, final int stepPosition, final int internalId) {
        hashCall(hashForward, internalId, stepPosition);
        hashCall(hashReverse, internalId + 1, stepPosition);
      }
    };
  }

  @Override
  public final void testReadLength() {
  }

  private static final class TestLongIncrementalHashLoop extends SearchIncrementalHashLoop {
    TestLongIncrementalHashLoop(int stepSize, final ExactHashFunction function) {
      super(stepSize, function, getOutputVarsEmpty(), null, null, false);
    }

    @Override
    public void hashCall(final long hash, final int internalId, final int stepPosition) {
    }
  }

  @Override
  protected void getLongLoop() throws IOException {
    final HashLoop hashLoop = new TestLongIncrementalHashLoop(1, new ExactHashFunction(1, 2));
    final ReaderParams re = new MockReaderParams(new ReaderLongMock(Integer.MAX_VALUE), SequenceMode.UNIDIRECTIONAL);
    final ISequenceParams se = new MockSequenceParams(re , 0, 0);
    hashLoop.execLoop(se, HashLoop.makeBuffer(se.reader()));
  }

  @Override
  public void testExceptionMessage() {
  }
}
