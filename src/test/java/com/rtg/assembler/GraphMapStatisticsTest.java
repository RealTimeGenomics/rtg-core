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

package com.rtg.assembler;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class GraphMapStatisticsTest extends TestCase {
  void inc(GraphMapStatistics stats, GraphMapStatistics.Stat s, int amount) {
    for (int i = 0; i < amount; ++i) {
      stats.increment(s);
    }
  }
  public void test() {
    final GraphMapStatistics stats = new GraphMapStatistics(null);
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
    final GraphMapStatistics combined = new GraphMapStatistics(null);
    final GraphMapStatistics first = new GraphMapStatistics(null);
    inc(first, GraphMapStatistics.Stat.AVOIDED_ALIGNMENTS, 3);
    inc(first, GraphMapStatistics.Stat.CROSS_CONTIG_SINGLE, 30);
    inc(first, GraphMapStatistics.Stat.SINGLE_END, 40);
    inc(first, GraphMapStatistics.Stat.PAIRED_END, 500);
    final GraphMapStatistics second = new GraphMapStatistics(null);
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
