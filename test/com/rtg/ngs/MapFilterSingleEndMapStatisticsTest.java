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
public class MapFilterSingleEndMapStatisticsTest extends TestCase {

  public void testIncrement() {
    final MapFilterSingleEndMapStatistics stats = new MapFilterSingleEndMapStatistics(null);
    stats.increment(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT);
    stats.set(MapStatisticsField.UNMATED_AMBIG_READS, Arm.RIGHT, 10);
    stats.increment(MapStatisticsField.UNMAPPED_UNMATED_POOR, Arm.RIGHT);
    stats.increment(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT);
    stats.set(MapStatisticsField.TOTAL_READS, Arm.LEFT, 15);
    stats.set(MapStatisticsField.TOTAL_READS, Arm.RIGHT, 15);
    final String str = stats.getStatistics();
    TestUtils.containsAll(str, "READ MAPPINGS",
                               "10  66.7% mapped",
                               " 1   6.7% unmapped with poor hits (XC = D)",
                               " 1   6.7% unmapped with no hits",
                               "15 100.0% total");

  }

  public void testBadFields() {
    final MapFilterSingleEndMapStatistics stats = new MapFilterSingleEndMapStatistics(null);
    try {
      stats.increment(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
    try {
      stats.set(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT, 10);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
    try {
      stats.value(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT);
      fail();
    } catch (final UnsupportedOperationException e) {
      //expected
    }
  }

}
