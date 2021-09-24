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
package com.rtg.index.hash.ngs.instances;

import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

public class DummyCG2MaskTest extends TestCase {

  private class MyAbstractCGMask extends AbstractCG2Mask {

    MyAbstractCGMask(final int windowSize, final ReadCall readCall, final TemplateCall templateCall) {
      super(windowSize, readCall, templateCall);

      assertEquals(5, cgAdjust3(5));

      final long o3 = 0b11111111110000000000000000L
                    | 0b00000001111111111111111111L;
      final long o4 = 0b1111111111000000000000000L
                    | 0b0000001111111111111111111L;
      assertEquals(0b11111111111111111111111111111L, cgAdjust3(o3));
      assertEquals(0b11111111111111111111111111111L, cgAdjust4(o4));
      assertEquals(0b01111111111111111111111111111L, cgAdjust3(o4));

      final long o3b = 0b11011011010000000000000000L
                     | 0b00000001011011011011011011L;
      final long o4b = 0b1101101101000000000000000L
                     | 0b0000001101101101101101101L;
      assertEquals(0b11011011011011011011011011011L, cgAdjust3(o3b));
      assertEquals(0b11011011011101101101101101101L, cgAdjust4(o4b));
    }

    @Override
    public void readAll(int readId, long v0, long v1) {
      // unused

    }
    @Override
    public void templateAll(int endPosition, long v0, long v1) {
      // unused

    }
  }

  public void testCgAdjust() {
    Diagnostic.setLogStream();
    new MyAbstractCGMask(18, null, null);
  }

}
