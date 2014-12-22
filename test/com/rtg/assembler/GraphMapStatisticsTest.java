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
public class GraphMapStatisticsTest extends TestCase {
  void inc(GraphMapStatistics stats, GraphMapStatistics.Stat s, int amount) {
    for (int i = 0; i < amount; i++) {
      stats.increment(s);
    }
  }
  public void test() {
    GraphMapStatistics stats = new GraphMapStatistics(null);
    inc(stats, GraphMapStatistics.Stat.AVOIDED_ALIGNMENTS, 3);
    inc(stats, GraphMapStatistics.Stat.CROSS_CONTIG_SINGLE, 30);
    inc(stats, GraphMapStatistics.Stat.SINGLE_END, 40);
    inc(stats, GraphMapStatistics.Stat.PAIRED_END, 500);
    TestUtils.containsAll(stats.getStatistics()
        , "  3 " + GraphMapStatistics.Stat.AVOIDED_ALIGNMENTS.mName
        , "500 " + GraphMapStatistics.Stat.PAIRED_END.mName
        , " 40 " + GraphMapStatistics.Stat.SINGLE_END.mName
        , " 30 " + GraphMapStatistics.Stat.CROSS_CONTIG_SINGLE.mName
    );
  }

  public void testCombine() {
    GraphMapStatistics combined = new GraphMapStatistics(null);
    GraphMapStatistics first = new GraphMapStatistics(null);
    inc(first, GraphMapStatistics.Stat.AVOIDED_ALIGNMENTS, 3);
    inc(first, GraphMapStatistics.Stat.CROSS_CONTIG_SINGLE, 30);
    inc(first, GraphMapStatistics.Stat.SINGLE_END, 40);
    inc(first, GraphMapStatistics.Stat.PAIRED_END, 500);
    GraphMapStatistics second = new GraphMapStatistics(null);
    inc(first, GraphMapStatistics.Stat.AVOIDED_ALIGNMENTS, 3);
    inc(first, GraphMapStatistics.Stat.CROSS_CONTIG_SINGLE, 30);
    inc(first, GraphMapStatistics.Stat.SINGLE_END, 40);
    inc(first, GraphMapStatistics.Stat.PAIRED_END, 500);
    combined.accumulate(first);
    combined.accumulate(new GraphMapStatistics(null));
    combined.accumulate(second);
    TestUtils.containsAll(combined.getStatistics()
        , "   6 " + GraphMapStatistics.Stat.AVOIDED_ALIGNMENTS.mName
        , "1000 " + GraphMapStatistics.Stat.PAIRED_END.mName
        , "  80 " + GraphMapStatistics.Stat.SINGLE_END.mName
        , "  60 " + GraphMapStatistics.Stat.CROSS_CONTIG_SINGLE.mName
    );

  }
}
