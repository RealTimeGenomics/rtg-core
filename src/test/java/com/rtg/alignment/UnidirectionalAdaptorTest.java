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
    assertArrayEquals(new int[] {1, 2, 3}, unidirectionalAdaptor.calculateEditDistance(new byte[0], 0, new byte[0], 0, false, 0, 0, false));
  }

}
