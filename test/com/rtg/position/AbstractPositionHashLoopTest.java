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
    final ReaderParams re = new MockReaderParams(sr, mode);
    final ISequenceParams se = new MockSequenceParams(re , 0, 2);
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
    final HashLoop loop = getHashLoop1a(15, 1, SequenceMode.UNIDIRECTIONAL.codeType().bits());
    final SequencesReader sr = new MockSequencesReader(SequenceType.DNA, Integer.MAX_VALUE + 1L);
    final ReaderParams re = new MockReaderParams(sr, SequenceMode.UNIDIRECTIONAL);
    final ISequenceParams se = new MockSequenceParams(re, 0, 1) {
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
