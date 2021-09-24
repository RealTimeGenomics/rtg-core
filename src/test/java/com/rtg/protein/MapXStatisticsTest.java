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
