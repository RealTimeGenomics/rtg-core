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

import com.rtg.reader.Arm;
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
