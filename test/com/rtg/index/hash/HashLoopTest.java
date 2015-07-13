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

    public MaxLenDummyReader(long maxSequenceLength) {
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

