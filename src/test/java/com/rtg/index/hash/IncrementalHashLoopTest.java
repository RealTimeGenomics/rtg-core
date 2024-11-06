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
package com.rtg.index.hash;

import java.io.IOException;

import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.SequencesReader;

/**
 */
public class IncrementalHashLoopTest extends AbstractIncrementalHashLoopTest {

  @Override
  protected HashLoop getHashLoop1(final int windowSize, int stepSize,
      final int maxId, final int bits, final int[] count) {
    final HashLoop hashLoop;
    final int[] mi = new int[1];
    mi[0] = maxId;
    hashLoop =
      new IncrementalHashLoop(stepSize, new ExactHashFunction(windowSize, bits), false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        assertTrue("" + internalId, internalId >= 0 && internalId < mi[0]);
        count[0]++;
      }

      @Override
      public void hashCallBidirectional(long hashForward, long hashReverse, int stepPosition, int internalId) {
        throw new UnsupportedOperationException("Not supported.");
      }
    };
    return hashLoop;
  }

  @Override
  protected HashLoop getHashLoop1a(final int windowSize, int stepSize, final int bits) {
    return new IncrementalHashLoop(stepSize, new ExactHashFunction(windowSize, bits), false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        // do nothing
      }

      @Override
      public void hashCallBidirectional(long hashForward, long hashReverse, int stepPosition, int internalId) {
        throw new UnsupportedOperationException("Not supported.");
      }
    };
  }

  @Override
  protected void checkCount(final SequencesReader sr, final int windowSize, final int stepSize, final SequenceMode mode, final int expected, final int maxId) throws IOException {
    final int bits = mode.codeType().bits();
    //System.err.println(" mode=" + mode + " type=" +  mode.type() + " bits=" + bits);
    final int[] count = new int[1];
    final HashLoop hashLoop = getHashLoop1(windowSize, stepSize, maxId, bits, count);
    final ReaderParams re = new MockReaderParams(sr);
    final ISequenceParams se = new MockSequenceParams(re, mode, 0, 1);
    hashLoop.execLoop(se, HashLoop.makeBuffer(sr));
    assertEquals(expected, count[0]);
    //test bad case when buffer supplied
    try {
      final HashLoop hashLoop1 = getHashLoop1a(windowSize, stepSize, bits);
      hashLoop1.execLoop(se, new byte[(int) Math.max(0, sr.maxLength() - 1)]);
    } catch (final RuntimeException e) {
      final String str = e.getMessage();
      assertTrue(str.startsWith("Allocated buffer too short. Allocated length="));
      assertTrue(str.contains(" Required length="));
    }
  }


  @Override
  protected HashLoop getHashLoop3(final int windowSize, int stepSize,
      final long[] expectedL, final int[] expectedI, final int bits, final int[] count) {
    final HashLoop hashLoop;
    hashLoop =
      new IncrementalHashLoop(stepSize, new ExactHashFunction(windowSize, bits), false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        assertEquals(count[0] + "", expectedL[count[0]], hash);
        assertEquals(count[0] + "", expectedI[count[0]], internalId);
        count[0]++;
      }

      @Override
      public void hashCallBidirectional(long hashForward, long hashReverse, int stepPosition, int internalId) {
        throw new UnsupportedOperationException("Not supported.");
      }
    };
    return hashLoop;
  }
}
