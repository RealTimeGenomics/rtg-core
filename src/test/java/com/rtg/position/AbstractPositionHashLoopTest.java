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

import java.io.IOException;

import com.rtg.index.hash.AbstractIncrementalHashLoopTest;
import com.rtg.index.hash.HashLoop;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.SequencesReader;
import com.rtg.util.StringUtils;
import com.rtg.util.array.ImmutableIntArray;

/**
 */
public abstract class AbstractPositionHashLoopTest extends AbstractIncrementalHashLoopTest {

  /**
   * Test method for {@link com.rtg.index.hash.IncrementalHashLoop}.
   * @throws IOException if an I/O error occurs.
   */
  public void testReadLength() throws IOException {
    final String str = ""
      + ">x1" + StringUtils.LS
      + "ACT" + StringUtils.LS
      + ">x2" + StringUtils.LS
      + "ACTG" + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
      checkCount4(sr, 1, 1, TRANSLATED, 6, 12);
      checkCount4(sr, 1, 1, UNIDIRECTIONAL, 7, 2);
      checkCount4(sr, 1, 1, BIDIRECTIONAL, 14, 4);

      checkCount4(sr, 2, 1, UNIDIRECTIONAL, 5, 2);
      checkCount4(sr, 2, 1, BIDIRECTIONAL, 10, 4);
      checkCount4(sr, 2, 1, TRANSLATED, 0, 12);

      checkCount4(sr, 2, 2, UNIDIRECTIONAL, 3, 2);
      checkCount4(sr, 2, 2, BIDIRECTIONAL, 6, 4);
      checkCount4(sr, 2, 2, TRANSLATED, 0, 12);

      checkCount4(sr, 5, 1, UNIDIRECTIONAL, 0, 2);
      checkCount4(sr, 5, 1, BIDIRECTIONAL, 0, 4);
      checkCount4(sr, 5, 1, TRANSLATED, 0, 12);
    }
  }

  protected void checkCount4(final SequencesReader sr, final int windowSize, final int stepSize, final SequenceMode mode, final int expected, final int maxId) throws IOException {
    final int bits = mode.codeType().bits();
    //System.err.println(" mode=" + mode + " type=" +  mode.type() + " bits=" + bits);
    final int[] count = new int[1];
    final HashLoop hashLoop = getHashLoop1(windowSize, stepSize, maxId, bits, count);
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
    final HashLoop loop = getHashLoop1a(15, 1, UNIDIRECTIONAL.codeType().bits());
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
