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

import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;
import com.rtg.variant.realign.AllPaths;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;
import com.rtg.variant.realign.ScoreFastUnderflow;
import com.rtg.variant.realign.ScoreFastUnderflowCG;

import junit.framework.TestCase;

/**
 */
public class ScoreInterfaceMemoTest extends TestCase {

  public void test() throws IOException {
    final ScoreInterfaceMemoInterface memo = new ScoreInterfaceMemo();
    final MachineErrorParamsBuilder builder = MachineErrorParams.builder();
    final RealignParams me = new RealignParamsImplementation(builder.create());
    final RealignParams me2 = new RealignParamsImplementation(builder.create());
    final AllPaths s1 = memo.getScoreInterface(me);
    assertTrue(s1 instanceof ScoreFastUnderflow);
    assertTrue(s1 == memo.getScoreInterface(me));

    final AllPaths s2 = memo.getScoreInterface(me2);
    assertTrue(s2 instanceof ScoreFastUnderflow);
    assertFalse(s1 == s2);

    final MachineErrorParams complete = builder.errors("complete").create();
    final RealignParams rpcg = new RealignParamsImplementation(complete);
    final AllPaths s3 = memo.getScoreInterface(rpcg);
    assertTrue(s3 instanceof ScoreFastUnderflowCG);
    assertTrue(s3 == memo.getScoreInterface(rpcg));
  }
}
