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

package com.rtg.assembler;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class PacBioStatisticsTest extends TestCase {
  public void testStatistics() {

  }
  void inc(PacBioStatistics stats, PacBioStatistics.Stat s, int amount) {
    for (int i = 0; i < amount; ++i) {
      stats.increment(s);
    }
  }
  public void test() {
    final PacBioStatistics stats = new PacBioStatistics(null);
    inc(stats, PacBioStatistics.Stat.TOTAL_READS, 3);
    inc(stats, PacBioStatistics.Stat.INTERNAL_READS, 500);
    inc(stats, PacBioStatistics.Stat.CROSS_CONTIG, 40);
    TestUtils.containsAll(stats.getStatistics()
        , "  3 " + PacBioStatistics.Stat.TOTAL_READS.mName
        , "500 " + PacBioStatistics.Stat.INTERNAL_READS.mName
        , " 40 " + PacBioStatistics.Stat.CROSS_CONTIG.mName
    );
  }
  public void testCombine() {
    final PacBioStatistics stats = new PacBioStatistics(null);
    inc(stats, PacBioStatistics.Stat.TOTAL_READS, 3);
    inc(stats, PacBioStatistics.Stat.INTERNAL_READS, 500);
    inc(stats, PacBioStatistics.Stat.CROSS_CONTIG, 40);
    final PacBioStatistics stats2 = new PacBioStatistics(null);
    inc(stats2, PacBioStatistics.Stat.TOTAL_READS, 6);
    inc(stats2, PacBioStatistics.Stat.INTERNAL_READS, 520);
    inc(stats2, PacBioStatistics.Stat.CROSS_CONTIG, 41);
    stats.accumulate(stats2);
    TestUtils.containsAll(stats.getStatistics()
        , "   9 " + PacBioStatistics.Stat.TOTAL_READS.mName
        , "1020 " + PacBioStatistics.Stat.INTERNAL_READS.mName
        , "  81 " + PacBioStatistics.Stat.CROSS_CONTIG.mName
    );

  }
}

