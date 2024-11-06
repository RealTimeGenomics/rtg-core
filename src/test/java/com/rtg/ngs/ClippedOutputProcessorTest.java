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
package com.rtg.ngs;


import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;

import junit.framework.TestCase;

/**
 */
public class ClippedOutputProcessorTest extends TestCase {

  long mTemplateId;
  long mTemplateStart;

 class MockOutputProcessor implements OutputProcessor {

    @Override
    public OutputProcessor threadClone(final HashingRegion clipRegion) {
      return new ClippedOutputProcessor(this, clipRegion);
    }

    @Override
    public void process(final long templateId, final String frame, final int readId, final int tStart, final int score, final int scoreIndel) {
      mTemplateId = templateId;
      mTemplateStart = tStart;
    }

    @Override
    public void close() {
    }
    @Override
    public void finish() {
    }

    @Override
    public void threadFinish() {
        close();
    }
  }

  public final void test() throws IOException {
    final OutputProcessor op = new MockOutputProcessor();
    final OutputProcessor op2 = op.threadClone(new HashingRegion(5, 10, 20, 5, -1, -1));

    op2.process(14, null, 1, 24, 0, 0); // Not clipped
    assertEquals(14, mTemplateId);
    assertEquals(24, mTemplateStart);

    op2.process(5, null, 1, 2, 0, 0);  // Should be clipped
    assertEquals(14, mTemplateId);
    assertEquals(24, mTemplateStart);

    op2.process(20, null, 1, 4, 0, 0);  // Not clipped
    assertEquals(20, mTemplateId);
    assertEquals(4, mTemplateStart);

  }

}
