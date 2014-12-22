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
package com.rtg.position.output;

import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.BuildParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.BidirectionalFrame;
import com.rtg.mode.SequenceMode;

import junit.framework.TestCase;

/**
 */
public class NgsWordOutputTest extends TestCase {

  public void test() throws IOException {
    final int[] value = new int[4];
    final OutputProcessor output = new OutputProcessor() {
      int mI;
      @Override
      public void process(long templateName, String frame, int readId, int tStart, int score, int scoreIndel) {
        value[mI++] = tStart;
      }

      @Override
      public void finish() {
      }

      @Override
      public OutputProcessor threadClone(HashingRegion region) {
        return null;
      }

      @Override
      public void threadFinish() {
      }

      @Override
      public void close() {
      }
    };
    final PositionParams params = PositionParams.builder()
            .buildParams(BuildParams.builder()
                    .windowSize(18).stepSize(18)
                    .sequences(SequenceParams.builder().mode(SequenceMode.UNIDIRECTIONAL).readerParam(new MockReaderParams(10, 10, SequenceMode.UNIDIRECTIONAL)).create())

                    .create())
            .searchParams(BuildParams.builder()
                    .windowSize(18).stepSize(1)
                    .sequences(SequenceParams.builder().mode(SequenceMode.BIDIRECTIONAL).readerParam(new MockReaderParams(10, 10, SequenceMode.BIDIRECTIONAL)).create())
                    .create())
            .create();
    final NgsWordOutput word = new NgsWordOutput(params, output);
    word.nextSequence(0, 10);
    word.nextQuery(BidirectionalFrame.FORWARD, 0);
    word.setPosition(10);
    word.hit(0, 1);
    word.hit(0, 2);
    word.hit(0, 3);
    word.hit(0, 4);
    assertEquals(10 - 1, value[0]);
    assertEquals(10 - 2, value[1]);
    assertEquals(10 - 3, value[2]);
    assertEquals(10 - 4, value[3]);
  }
}
