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
package com.rtg.protein;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.reader.Arm;
import com.rtg.ngs.MapStatisticsField;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;


/**
 */
public class MapXStatisticsTest extends TestCase {

  public void test() throws IOException {
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    final PrintStream ps = new PrintStream(log);
    Diagnostic.setLogStream(ps);
    final MapXStatistics stats = new MapXStatistics(null);
    try {
      stats.printStatistics(null);
    } finally {
      //closes log
      Diagnostic.setLogStream();
    }
    assertEquals("", log.toString());
    stats.set(MapStatisticsField.TOTAL_READS, Arm.LEFT, 100);
    stats.set(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT, 50);
    stats.set(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.LEFT, 20);
    stats.set(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT, 30);
    final ByteArrayOutputStream out = new ByteArrayOutputStream();
    stats.printStatistics(out);
    final String outString = out.toString();
    TestUtils.containsAll(outString,
        "ALIGNMENT STATISTICS",
        "100 100.0% total",
        " 50  50.0% alignments met threshold",
        " 20  20.0% alignments failed threshold",
        " 30  30.0% no significant alignments"
        );
  }
}
