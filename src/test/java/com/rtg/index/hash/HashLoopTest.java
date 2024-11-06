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

import com.rtg.launcher.ISequenceParams;
import com.rtg.mode.Frame;
import com.rtg.reader.DummySequencesReader;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class HashLoopTest extends TestCase {

  private class MyHashLoop extends HashLoop {

    /**
     * @param stepSize step size
     * @param windowSize windows size
     */
    protected MyHashLoop(int windowSize, int stepSize) {
      super(windowSize, stepSize);
      // do nothing
    }

    @Override
    public void end() {
      // do nothing
    }

    @Override
    public void endAll() {
      // do nothing
    }

    @Override
    public void endSequence() {
      // do nothing
    }

    @Override
    public long execLoop(ISequenceParams params, byte[] byteBuffer) {
      // do nothing
      return 0;
    }

    @Override
    public void hashCall(long hash, int internalId, int stepPosition) {
      // do nothing
    }

    @Override
    public void hashCallBidirectional(long hashForward, long hashReverse, int stepPosition, int internalId) {
      // do nothing
    }

    @Override
    public void next(long seq, Frame frame) {
      // do nothing
    }
  }

  public final void test() {
    final MyHashLoop hl = new MyHashLoop(20, 20);
    hl.setThreadPadding(20);
    assertEquals(20, hl.getThreadPadding());
  }

  private class MaxLenDummyReader extends DummySequencesReader {

    private final long mMaxSequenceLength;

    MaxLenDummyReader(long maxSequenceLength) {
      mMaxSequenceLength = maxSequenceLength;
    }

    @Override
    public long maxLength() {
      return mMaxSequenceLength;
    }
  }

  public final void testMakeBuffer() {
    try (MemoryPrintStream stream = new MemoryPrintStream()) {
      try {
        Diagnostic.setLogStream(stream.printStream());
        try {
          HashLoop.makeBuffer(new MaxLenDummyReader(Integer.MAX_VALUE));
          fail();
        } catch (final SlimException e) {
          //expected
        }
        assertTrue(stream.toString().contains("There is a sequence which is too long to process. Its length is \"" + Integer.MAX_VALUE + "\" bytes. See the SDF output for the name of the sequence."));
      } finally {
        Diagnostic.setLogStream();
      }
      final byte[] buffer = HashLoop.makeBuffer(new MaxLenDummyReader(3));
      assertEquals(3, buffer.length);
    }
  }

}

