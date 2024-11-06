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

