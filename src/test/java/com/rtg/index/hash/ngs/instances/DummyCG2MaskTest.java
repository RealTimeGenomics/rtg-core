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
