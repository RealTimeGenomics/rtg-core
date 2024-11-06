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
    final ReaderParams re = new MockReaderParams(new ReaderLongMock(Integer.MAX_VALUE));
    final ISequenceParams se = new MockSequenceParams(re, SequenceMode.UNIDIRECTIONAL, 0, 0);
    hashLoop.execLoop(se, HashLoop.makeBuffer(se.reader()));
  }

  @Override
  public void testExceptionMessage() {
  }
}
