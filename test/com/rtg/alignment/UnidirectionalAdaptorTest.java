/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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

import static org.junit.Assert.assertArrayEquals;

import org.junit.Test;

/**
 */
public class UnidirectionalAdaptorTest {
  private static class MockUnidirectionalEditDistance implements UnidirectionalEditDistance {
    @Override
    public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
      return new int[] {1, 2, 3};
    }

    @Override
    public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int maxScore, int maxShift) {
      throw new UnsupportedOperationException();
    }

    @Override
    public int[] calculateEditDistanceFixedEnd(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateExpectedStartPos, int templateEndPos, int maxScore, int maxShift) {
      throw new UnsupportedOperationException();
    }

    @Override
    public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int templateEndPos, int maxScore, int maxShift) {
      throw new UnsupportedOperationException();
    }

    @Override
    public void logStats() {
      throw new UnsupportedOperationException();
    }

  }

  @Test
  public void adaptorTest() {
    final UnidirectionalAdaptor unidirectionalAdaptor = new UnidirectionalAdaptor(new MockUnidirectionalEditDistance());
    assertArrayEquals( new int[] {1, 2, 3}, unidirectionalAdaptor.calculateEditDistance(new byte[0], 0, new byte[0], 0, false, 0, 0, false));
  }

}