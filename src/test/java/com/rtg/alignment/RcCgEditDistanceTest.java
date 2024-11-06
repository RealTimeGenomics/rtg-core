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

import java.io.IOException;

import com.rtg.AbstractTest;
import com.rtg.mode.DnaUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.variant.MachineErrorParamsBuilder;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;

/**
 */
public class RcCgEditDistanceTest extends AbstractTest {

  public void testAbstract() {
    final CgGotohEditDistanceTest test = new CgGotohEditDistanceTest();
    test.testCGOverlap0Rev();
    test.testCgOverlapRealWorld7();
  }

  private static class TestGotohEditDistance extends CgGotohEditDistance {
    TestGotohEditDistance(int maxShift, RealignParams params) {
      super(maxShift, params, 0);
    }
    @Override
    public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
      return ActionsHelper.build("=====B====X===============NNNNNN==========", 0, 2);
    }
  }

  public void testCGGotohHandling() {
    final RealignParamsImplementation mParams;
    try {
      mParams = new RealignParamsImplementation(new MachineErrorParamsBuilder().errors("cg_test_errors").create());
    } catch (final InvalidParamsException | IOException e) {
      throw new RuntimeException(e);
    }

    final BidirectionalEditDistance rced = new RcCgEditDistance(new TestGotohEditDistance(7, mParams));
    final int[] actions = rced.calculateEditDistance(DnaUtils.encodeString("attcttactanccccaactgagccc     agtaggagta".replaceAll(" ", "")), 35, DnaUtils.encodeString("atactcctacttttgctgggctcagttgggggaagtagaat"), 0, true, 25, 7, true);

    assertEquals(1, ActionsHelper.zeroBasedTemplateStart(actions));

  }
}
