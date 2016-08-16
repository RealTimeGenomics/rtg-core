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
        "ED#         Null#   MaxOk   MaxBad  TooLow  Correct TooHigh Total   AvUsecs".replaceAll("  *", "\t"),
        /*NoIndel*/ " 41      0       0       0       0      0       41 ".replaceAll("  *", "\t"),
        /*LowerB*/ " 41      0      0       0       0       0       41 ".replaceAll("  *", "\t"),
        /*HopSte*/ " 41      0       0       0       0      0       41 ".replaceAll("  *", "\t"),
        /*Seeded*/ " 0      0       0       0       41      0       41 ".replaceAll("  *", "\t"),
        /*GotohE*/ " 0       0      0       0       41      0       41 ".replaceAll("  *", "\t")
        );
  }
}
