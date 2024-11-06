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

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;

import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.util.Resources;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * This tests the whole alignment chain.
 */
public class EditDistanceChainTest extends TestCase {

  private EditDistanceChainAnalyzer mEd = null;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(1).create();
    mEd = new EditDistanceChainAnalyzer(params);
  }

  @Override
  public void tearDown() {
    mEd = null;
  }

  public void testGetName() {
    int i = 0;
    assertEquals(NoIndelsEditDistance.class.getSimpleName(), mEd.getName(i++));
    assertEquals(SingleIndelEditDistance.class.getSimpleName(), mEd.getName(i++));
    assertEquals(SingleIndelSeededEditDistance.class.getSimpleName(), mEd.getName(i++));
    assertEquals(LowerBoundEditDistance.class.getSimpleName(), mEd.getName(i++));
    assertEquals(HopStepEditDistanceLong.class.getSimpleName(), mEd.getName(i));
  }

  public void test1000a() throws IOException {
    final InputStream input = Resources.getResourceAsStream("com/rtg/alignment/resources/read1000a.txt");
    assertNotNull(input);
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    Diagnostic.setLogStream(new PrintStream(bos));
    mEd.setShift(20);
    //mEd.printEachRead(true);
    //mEd.printEachError(true);
    mEd.doFile(new BufferedReader(new InputStreamReader(input)), 50);
    mEd.logStats();
    Diagnostic.setLogStream();
    final String stats = bos.toString();
//    System.err.println(stats);
    TestUtils.containsAll(stats,
        "EditDistanceChainTest statistics, readlen: 1000 maxScore: 50",
        "ED#         Null#   MaxOk   MaxBad  TooLow  Correct TooHigh Total   AvUsecs".replaceAll(" +", "\t"),
        /*NoIndel*/ " 41      0       0       0       0      0       41 ".replaceAll(" +", "\t"),
        /*LowerB*/ " 41      0      0       0       0       0       41 ".replaceAll(" +", "\t"),
        /*HopSte*/ " 41      0       0       0       0      0       41 ".replaceAll(" +", "\t"),
        /*Seeded*/ " 0      0       0       0       41      0       41 ".replaceAll(" +", "\t"),
        /*GotohE*/ " 0       0      0       0       41      0       41 ".replaceAll(" +", "\t")
        );
  }
}
