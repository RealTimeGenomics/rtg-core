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

import com.rtg.index.Index;
import com.rtg.index.hash.ExactHashFunction;
import com.rtg.index.hash.HashFunction;
import com.rtg.index.hash.HashLoop;
import com.rtg.index.hash.ResetHashLoop;
import com.rtg.index.hash.ResetHashLoopTest;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.SequencesReader;
import com.rtg.util.StringUtils;

/**
 */
public class SearchResetHashLoopTest extends ResetHashLoopTest {

  private static FinderPositionOutput getOutputVarsEmpty() {
    return new FinderPositionOutput(null, new MockPositionOutput());
  }

  private static final class TestLoopBuildResetHashLoop extends SearchResetHashLoop {
    TestLoopBuildResetHashLoop(final int windowSize, int stepSize, final HashFunction function) {
      super(windowSize, stepSize, function, getOutputVarsEmpty()/*output*/, null/*finder*/, null/*index*/, false);
    }

    @Override
    public void hashCall(final long hash, final int internalId, final int stepPosition) {
    }
  }

  @Override
  protected HashLoop getHashLoop(final int windowSize, int stepSize, final HashFunction function) throws IOException {
    return new TestLoopBuildResetHashLoop(windowSize, stepSize, function);
  }

  @Override
  protected ResetHashLoop getHashLoop1a(final int windowSize, int stepSize, final int bits)
  throws IOException {
    return new SearchResetHashLoop(windowSize,  stepSize, new ExactHashFunction(windowSize, bits), getOutputVarsEmpty(), (Index) null, false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) {
        // do nothing
      }
    };
  }

  @Override
  protected HashLoop getHashLoop1(final int windowSize, int stepSize,
      final int maxId, final int bits, final int[] count) throws IOException {
    final HashLoop hashLoop;
    final int[] mi = new int[1];
    mi[0] = maxId;
    hashLoop =
      new SearchResetHashLoop(windowSize,  stepSize, new ExactHashFunction(windowSize, bits), getOutputVarsEmpty(), (Index) null, false) {
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
      final long[] expectedL, final int[] expectedI, final int bits) throws IOException {
    final HashLoop hashLoop;
    final int[] count = new int[1];
    //System.err.println(" mode=" + mode + " type=" +  mode.type() + " bits=" + bits);
    hashLoop =
      new SearchResetHashLoop(windowSize,  stepSize, new ExactHashFunction(windowSize, bits), getOutputVarsEmpty(), (Index) null, false) {
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
      + ">x1" + StringUtils.LS
      + "ACT" + StringUtils.LS
      + ">x2" + StringUtils.LS
      + "ACTG" + StringUtils.LS;
    try (SequencesReader sr = getReaderDNA(str)) {
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
    final Index index = new MockIndex();
    final int[] count = new int[1];
    final int[] mi = new int[1];
    mi[0] = maxId;
    final HashLoop hashLoop =
      new SearchResetHashLoop(windowSize,  stepSize, new ExactHashFunction(windowSize, bits), getOutputVarsEmpty(), index, false) {
      @Override
      public void hashCall(final long hash, final int internalId, final int stepPosition) throws IOException {
        assertTrue(" id=" + internalId + " max=" + mi[0], internalId >= 0 && internalId < mi[0]);
        count[0]++;
        //System.err.println(" internalId=" + internalId + " stepPosition=" + stepPosition);
        super.hashCall(hash, internalId, stepPosition);
      }
    };
    final ReaderParams re = new MockReaderParams(sr, mode);
    final ISequenceParams se = new MockSequenceParams(re , 0, 2);
    hashLoop.execLoop(se, HashLoop.makeBuffer(sr));
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
}
