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

import static com.rtg.mode.SequenceMode.BIDIRECTIONAL;
import static com.rtg.mode.SequenceMode.TRANSLATED;
import static com.rtg.mode.SequenceMode.UNIDIRECTIONAL;
import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import com.rtg.index.Index;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.index.hash.HashFunction;
import com.rtg.index.hash.HashLoop;
import com.rtg.index.hash.ResetHashLoop;
import com.rtg.index.hash.ResetHashLoopTest;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.SequencesReader;
import com.rtg.util.array.ImmutableIntArray;

/**
 */
public class BuildResetHashLoopTest extends ResetHashLoopTest {

  final class TestLoopBuildResetHashLoop extends BuildResetHashLoop {
    TestLoopBuildResetHashLoop(final int windowSize, int stepSize, final HashFunction function) {
      super(windowSize, stepSize, function, null, 0/*mxs*/);
    }

    @Override
    public void hashCall(final long hash, final int internalId, final int stepPosition) {
    }
  }

  @Override
  protected HashLoop getHashLoop(final int windowSize, int stepSize, final HashFunction function) {
    return new TestLoopBuildResetHashLoop(windowSize, stepSize, function);
  }

  @Override
  protected ResetHashLoop getHashLoop1a(final int windowSize, int stepSize, final int bits) {
    return new BuildResetHashLoop(windowSize,  stepSize, new ExactHashFunction(windowSize, bits), null/*index*/, 0/*mxs*/) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        // do nothing
      }
      @Override
      public void hashCallBidirectional(long hashForward, long hashReverse, int stepPosition, int internalId) {
      }
    };
  }

  @Override
  protected HashLoop getHashLoop1(final int windowSize, int stepSize,
      final int maxId, final int bits, final int[] count) {
    final HashLoop hashLoop;
    final int[] mi = new int[1];
    mi[0] = maxId;
    hashLoop =
      new BuildResetHashLoop(windowSize,  stepSize, new ExactHashFunction(windowSize, bits), null, 0/*mxs*/) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        assertTrue("" + internalId, internalId >= 0 && internalId < mi[0]);
        count[0]++;
      }
    };
    return hashLoop;
  }

  @Override
  protected HashLoop getHashLoop3(final int windowSize, int stepSize,
      final long[] expectedL, final int[] expectedI, final int bits) {
    final HashLoop hashLoop;
    final int[] count = new int[1];
    //System.err.println(" mode=" + mode + " type=" +  mode.type() + " bits=" + bits);
    hashLoop =
      new BuildResetHashLoop(windowSize,  stepSize, new ExactHashFunction(windowSize, bits), null, 0/*mxs*/) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        assertEquals(count[0] + "", expectedL[count[0]], hash);
        assertEquals(count[0] + "", expectedI[count[0]], internalId);
        count[0]++;
      }
    };
    return hashLoop;
  }

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public final void testReadLength() throws IOException {
    final String str = ""
      + ">x1" + LS
      + "ACT" + LS
      + ">x2" + LS
      + "ACTG" + LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      final String exp1 = ""
        + "0" + LS
        + "126" + LS
        + "252" + LS
        + "294" + LS
        + "378" + LS
        + "420" + LS;
      checkCount4(sr, 1, 1, TRANSLATED, 6, 12, exp1);
      final String exp2 = ""
        + "0" + LS
        + "1" + LS
        + "2" + LS
        + "42" + LS
        + "43" + LS
        + "44" + LS
        + "45" + LS;
      checkCount4(sr, 1, 1, UNIDIRECTIONAL, 7, 2, exp2);
      checkCount4(sr, 1, 1, BIDIRECTIONAL, 14, 4, null);

      checkCount4(sr, 2, 1, UNIDIRECTIONAL, 5, 2, null);
      checkCount4(sr, 2, 1, BIDIRECTIONAL, 10, 4, null);
      checkCount4(sr, 2, 1, TRANSLATED, 0, 12, null);

      checkCount4(sr, 2, 2, UNIDIRECTIONAL, 3, 2, null);
      checkCount4(sr, 2, 2, BIDIRECTIONAL, 6, 4, null);
      checkCount4(sr, 2, 2, TRANSLATED, 0, 12, null);

      checkCount4(sr, 5, 1, UNIDIRECTIONAL, 0, 2, null);
      checkCount4(sr, 5, 1, BIDIRECTIONAL, 0, 4, null);
      checkCount4(sr, 5, 1, TRANSLATED, 0, 12, null);
    }
  }

  protected void checkCount4(final SequencesReader sr, final int windowSize, final int stepSize, final SequenceMode mode, final int expected, final int maxId, final String expId) throws IOException {
    final int bits = mode.codeType().bits();
    //System.err.println(" mode=" + mode + " type=" +  mode.type() + " bits=" + bits);
    final Index index = new MockIndex();
    final int[] count = new int[1];
    final int[] mi = new int[1];
    mi[0] = maxId;
    final HashLoop hashLoop =
      new BuildResetHashLoop(windowSize,  stepSize, new ExactHashFunction(windowSize, bits), index, 42) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        assertTrue(" id=" + internalId + " max=" + mi[0], internalId >= 0 && internalId < mi[0]);
        count[0]++;
        //System.err.println(" internalId=" + internalId + " stepPosition=" + stepPosition);
        super.hashCall(hash, internalId, stepPosition);
      }
    };
    final ReaderParams re = new MockReaderParams(sr);
    final ISequenceParams se = new MockSequenceParams(re, mode, 0, 2);
    try {
      hashLoop.readLengths();
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Read lengths not yet constructed", e.getMessage());
    }
    hashLoop.execLoop(se, HashLoop.makeBuffer(sr));
    final ImmutableIntArray rl = hashLoop.readLengths();
    assertEquals("[3, 4]", rl.toString());
    assertEquals(expected, count[0]);
    //System.err.println(index.toString());
    if (expId != null) {
      assertEquals(expId, index.toString());
    }
  }


  /**
   * Test method for {@link BuildIncrementalHashLoop#nextSeq(int, int)}.
   */
  public final void testNextSeq() {

  }

  /**
   * Test method for {@link BuildIncrementalHashLoop#readLengths()}.
   */
  public final void testReadLengths() {
  }

  public void testExceptionMessage() throws Exception {
    final HashLoop loop = new BuildResetHashLoop(15, 1, new ExactHashFunction(15, UNIDIRECTIONAL.codeType().bits()), null, 42) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        // do nothing
      }
    };
    final SequencesReader sr = new MockSequencesReader(SequenceType.DNA, Integer.MAX_VALUE + 1L);
    final ReaderParams re = new MockReaderParams(sr);
    final ISequenceParams se = new MockSequenceParams(re, UNIDIRECTIONAL, 0, 1) {
      @Override
      public HashingRegion region() {
        return new HashingRegion(0, Integer.MAX_VALUE + 1L);
      }
    };
    try {
      loop.execLoop(se, null);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("too many reads:2147483648", e.getMessage());
    }
  }

}
