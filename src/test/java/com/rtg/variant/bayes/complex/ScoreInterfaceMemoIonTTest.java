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
package com.rtg.variant.bayes.complex;

import java.io.IOException;
import java.io.StringReader;

import com.rtg.util.InvalidParamsException;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.realign.HomoPolymerParams;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ScoreInterfaceMemoIonTTest extends TestCase {

  public void test() throws IOException, InvalidParamsException {
    final ScoreInterfaceMemoIonT memo = new ScoreInterfaceMemoIonT(new HomoPolymerParams(SimplePossibility.SINGLETON, 2, 2, new StringReader("")), new HomoPolymerParams(LogApproximatePossibility.SINGLETON, 2, 2, new StringReader("")));
    final MachineErrorParams me = MachineErrorParams.builder().create();
    final RealignParams rp = new RealignParamsImplementation(me);
    assertTrue(memo.getScoreInterface(rp) == memo.getScoreInterface(rp));
    try {
      memo.getScoreInterface(new RealignParamsImplementation(MachineErrorParams.builder().errors("complete").create()));
      fail();
    } catch (final UnsupportedOperationException e) {

    }
  }
}
