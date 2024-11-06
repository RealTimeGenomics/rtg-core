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

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import com.rtg.reader.Arm;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;


/**
 */
public class SingleEndMapStatisticsTest extends TestCase {

  public void testPrintStatistics() throws Exception {
    final SingleEndMapStatistics testStats = new SingleEndMapStatistics(null);
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(log);
    Diagnostic.setLogStream(ps);
    try {
      testStats.printStatistics(null);
    } finally {
      //closes log
      Diagnostic.setLogStream();
    }
    assertEquals("", log.toString());
    testStats.set(MapStatisticsField.TOTAL_READS, Arm.LEFT, 1234L);
    testStats.set(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT, 567L);
    testStats.set(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT, 5L);
    testStats.set(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT, 932L);
    testStats.set(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT, 152L);
    testStats.set(MapStatisticsField.MISSING, Arm.LEFT, 0L);

    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    testStats.printStatistics(out);
    final String outString = out.toString();
    TestUtils.containsAll(outString,
        "READ MAPPINGS",
        "1234 100.0% total",
        " 932  75.5% mapped uniquely",
        "   5   0.4% mapped ambiguously",
        "   0   0.0% unmapped with too many hits (XC = C)",
        " 152  12.3% unmapped with no hits"
        );
    assertFalse(outString.contains("left arms missing"));
    assertFalse(outString.contains("right arms missing"));
  }

  public void testMisc() {
    final SingleEndMapStatistics testStats = new SingleEndMapStatistics(null);

    testStats.increment(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT);
    testStats.increment(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT);
    testStats.increment(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT);
    testStats.increment(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT);
    testStats.increment(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT);
    testStats.increment(MapStatisticsField.UNMAPPED_BLOCKED, Arm.LEFT);
    testStats.increment(MapStatisticsField.UNMAPPED_MATED_POOR, Arm.LEFT);
    testStats.increment(MapStatisticsField.UNMAPPED_MATED_TOO_MANY, Arm.LEFT);
    testStats.increment(MapStatisticsField.UNMAPPED_TOPN, Arm.LEFT);
    testStats.increment(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.LEFT);
    testStats.increment(MapStatisticsField.UNMAPPED_UNMATED_TOO_MANY, Arm.LEFT);
    testStats.increment(MapStatisticsField.MISSING, Arm.LEFT);
    testStats.increment(MapStatisticsField.TOTAL_READS, Arm.LEFT);

    assertEquals(1, testStats.value(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.UNMAPPED_BLOCKED, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.UNMAPPED_MATED_POOR, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.UNMAPPED_MATED_TOO_MANY, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.UNMAPPED_TOPN, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.UNMAPPED_UNMATED_TOO_MANY, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.MISSING, Arm.LEFT));
    assertEquals(1, testStats.value(MapStatisticsField.TOTAL_READS, Arm.LEFT));

    TestUtils.containsAll(testStats.getStatistics(),
        "READ MAPPINGS",
        "1 100.0% mapped uniquely (NH = 1)",
        "1 100.0% mapped ambiguously (NH > 1)",
        "1 100.0% unmapped due to read frequency (XC = B)",
        "1 100.0% unmapped with too many hits (XC = C)",
        "1 100.0% unmapped with poor hits (XC = D)",
        "1 100.0% unmapped with too many good hits (XC = E)",
        "1 100.0% unmapped with no hits",
        "1 100.0% total"
      );

    testStats.reset();

    assertEquals(0, testStats.value(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.UNMAPPED_BLOCKED, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.UNMAPPED_MATED_POOR, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.UNMAPPED_MATED_TOO_MANY, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.UNMAPPED_TOPN, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.UNMAPPED_UNMATED_TOO_MANY, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.MISSING, Arm.LEFT));
    assertEquals(0, testStats.value(MapStatisticsField.TOTAL_READS, Arm.LEFT));
  }
}
