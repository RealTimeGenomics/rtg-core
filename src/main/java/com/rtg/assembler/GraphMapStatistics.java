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

import java.io.File;
import java.util.EnumMap;
import java.util.Map;

import com.rtg.launcher.AbstractStatistics;
import com.rtg.util.StringUtils;

/**
 */
public class GraphMapStatistics extends AbstractStatistics {

  static class Counter {
    long mCount = 0;
    public void increment() {
      ++mCount;
    }
    public long getCount() {
      return mCount;
    }
  }
  /** Statistics represented in this object */
  public enum Stat {
    /** Number of single end reads */
    SINGLE_END("Single end reads"),
    /** Number of paired end reads */
    PAIRED_END("Paired end reads"),
    /** There were too many possible paths for this alignment */
    TOO_MANY_PATHS("Too many paths"),
    /** There were too many possible paths for this pair */
    TOO_MANY_PAIR_PATHS("Too many pairings"),
    /** There were no paths for this pair */
    NO_PAIRINGS("No pairings found"),
    /** No valid paths could be found */
    NO_PATHS("No paths"),
    /** Reads that were mapped individually */
    MAPPED("Mapped without pairing"),
    /** Reads which had no hits */
    NO_HITS("No hits"),
    /** Successful pairings */
    PAIRED("Successfully paired"),
    /** Alignments avoided due to equivalent hits */
    AVOIDED_ALIGNMENTS("Equivalent alignments skipped"),
    /** Paired alignments within a single contig */
    SINGLE_CONTIG_PAIR("Paired in a single contig"),
    /** Unpaired alignments within a single contig */
    SINGLE_CONTIG_MAPPING("Mapped in a single contig"),
    /** Paired alignments that cross contig boundaries */
    CROSS_CONTIG_PAIR("Paired across contigs"),
    /** Unpaired alignments that cross contig boundaries */
    CROSS_CONTIG_SINGLE("Reads mapped across contigs");

    final String mName;
    Stat(String name) {
      mName = name;
    }
  }

  final Map<Stat, Counter> mMap = new EnumMap<>(Stat.class);

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   */
  GraphMapStatistics(File outputDirectory) {
    super(outputDirectory);
    for (Stat s : Stat.values()) {
      mMap.put(s, new Counter());
    }
  }

  void increment(Stat s) {
    mMap.get(s).increment();
  }

  @Override
  protected String getStatistics() {
    int formatLength = 0;
    for (Stat s : Stat.values()) {
      formatLength = Math.max(formatLength, String.valueOf(mMap.get(s).getCount()).length());
    }
    final StringBuilder sb = new StringBuilder();
    for (Stat s : Stat.values()) {
      sb.append(String.format("%" + formatLength + "d %s", mMap.get(s).getCount(), s.mName));
      sb.append(StringUtils.LS);
    }
    return sb.toString();
  }

  /**
   * Add the contained statistics into this one
   * @param stats a set of statistics to add
   */
  public void accumulate(GraphMapStatistics stats) {
    for (Map.Entry<Stat, Counter> entry : stats.mMap.entrySet()) {
      mMap.get(entry.getKey()).mCount += entry.getValue().mCount;
    }
  }

  @Override
  public void generateReport() {
  }
}
