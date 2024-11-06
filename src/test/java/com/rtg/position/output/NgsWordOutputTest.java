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
                    .sequences(SequenceParams.builder().mode(SequenceMode.UNIDIRECTIONAL).readerParam(new MockReaderParams(10, 10, SequenceMode.UNIDIRECTIONAL.codeType())).create())

                    .create())
            .searchParams(BuildParams.builder()
                    .windowSize(18).stepSize(1)
                    .sequences(SequenceParams.builder().mode(SequenceMode.BIDIRECTIONAL).readerParam(new MockReaderParams(10, 10, SequenceMode.BIDIRECTIONAL.codeType())).create())
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
