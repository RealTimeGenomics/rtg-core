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
package com.rtg.ngs;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

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

  public void testMerge() {
    Diagnostic.setLogStream();
    final SingleEndMapStatistics testStats = new SingleEndMapStatistics(null);
    testStats.set(MapStatisticsField.TOTAL_READS, Arm.LEFT, 220L);
    testStats.set(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT, 25L);
    testStats.set(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT, 75L);
    testStats.set(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT, 80L);
    testStats.set(MapStatisticsField.MISSING, Arm.LEFT, 0L);

    final SingleEndMapStatistics testStats2 = new SingleEndMapStatistics(null);
    testStats2.set(MapStatisticsField.TOTAL_READS, Arm.LEFT, 220L);
    testStats2.set(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT, 80L);
    testStats2.set(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT, 20L);
    testStats2.set(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT, 7L);
    testStats2.set(MapStatisticsField.MISSING, Arm.LEFT, 2L);

    testStats.merge(testStats2);

    assertEquals(220L, testStats.value(MapStatisticsField.TOTAL_READS, Arm.LEFT));
    assertEquals(105L, testStats.value(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT));
    assertEquals(95L, testStats.value(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT));
    assertEquals(20L, testStats.value(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT));
    assertEquals(2L, testStats.value(MapStatisticsField.MISSING, Arm.LEFT));
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
