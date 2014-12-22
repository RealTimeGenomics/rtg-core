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
      mCount++;
    }
    public long getCount() {
      return mCount;
    }
  }
  /** Statistics represented in this object */
  public static enum Stat {
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
      formatLength = Math.max(formatLength, (mMap.get(s).getCount() + "").length());
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
