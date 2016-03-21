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

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 * Test class
 */
public class MapFilterPairedMapStatisticsTest extends TestCase {

  public void testIncrement() {
    final MapFilterPairedMapStatistics stats = new MapFilterPairedMapStatistics(null);
    stats.increment(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT);
    stats.set(MapStatisticsField.UNMATED_AMBIG_READS, Arm.RIGHT, 10);
    stats.increment(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.RIGHT);
    stats.increment(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT);
    stats.incrementBothUnmapped();
    stats.set(MapStatisticsField.TOTAL_READS, Arm.LEFT, 15);
    stats.set(MapStatisticsField.TOTAL_READS, Arm.RIGHT, 15);
    final String str = stats.getStatistics(true);
    TestUtils.containsAll(str, "ARM MAPPINGS",
                          " left right  both",
                          "    1    10    11  36.7% mapped",
                          "    0     1     1   3.3% unmapped with poor hits (XC = D)",
                          "    1     0     1   3.3% unmapped with no hits",
                          "   15    15    30 100.0% total");
  }

  public void testBadFields() {
    final MapFilterPairedMapStatistics stats = new MapFilterPairedMapStatistics(null);
    try {
      stats.increment(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
    try {
      stats.set(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT, 10);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
    try {
      stats.value(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
    try {
      stats.valueAsPercent(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
  }

}
