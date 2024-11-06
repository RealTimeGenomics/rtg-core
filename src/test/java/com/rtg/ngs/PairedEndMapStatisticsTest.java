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

import com.rtg.AbstractTest;
import com.rtg.reader.Arm;
import com.rtg.util.TestUtils;

/**
 */
public class PairedEndMapStatisticsTest extends AbstractTest {

  public void testPrintStatistics() throws Exception {
    final PairedEndMapStatistics testStats = new PairedEndMapStatistics(true, null);
    testStats.printStatistics(null);

    testStats.set(MapStatisticsField.TOTAL_READS, Arm.LEFT, 1000L);
    testStats.set(MapStatisticsField.TOTAL_READS, Arm.RIGHT, 234L);
    testStats.set(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.RIGHT, 567L);
    testStats.set(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT, 5L);
    testStats.set(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.RIGHT, 932L);
    testStats.set(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT, 152L);
    testStats.set(MapStatisticsField.MISSING, Arm.LEFT, 12L);
    testStats.set(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT, 33L);
    testStats.set(MapStatisticsField.MATED_AMBIG_READS, Arm.RIGHT, 66L);

    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    testStats.printStatistics(out);
    final String outString = out.toString();
    TestUtils.containsAll(outString,
        "ARM MAPPINGS",
        " left right  both",
        " 1000   234  1234 100.0% total",
        "    0     0     0   0.0% mated uniquely",
        "   33    66    99   8.0% mated ambiguously (NH > 1)",
        "    0   932   932  75.5% unmated uniquely",
        "    5     0     5   0.4% unmated ambiguously",
        "    0     0     0   0.0% unmapped with no matings but too many hits (XC = C)",
        "  152     0   152  12.3% unmapped with no hits",
        "    0     0     0   0.0% unmapped due to read frequency",
        "    0     0     0   0.0% unmapped with no matings",
        "    0     0     0   0.0% unmapped with poor matings",
        "    0     0     0   0.0% unmapped with too many matings",
        "    0     0     0   0.0% unmapped with no matings and poor hits",
        "    0     0     0   0.0% unmapped with no matings and too many good hits"
        );
  }
}
