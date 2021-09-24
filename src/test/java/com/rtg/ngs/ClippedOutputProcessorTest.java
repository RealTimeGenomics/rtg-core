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
