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

import static com.rtg.assembler.GraphMapStatistics.Counter;

import java.io.File;
import java.util.EnumMap;
import java.util.Map;

import com.rtg.launcher.AbstractStatistics;
import com.rtg.util.StringUtils;

/**
 */
public class PacBioStatistics extends AbstractStatistics  {
  /** Statistics represented in this object */
  public enum Stat {
    /** total number of reads */
    TOTAL_READS("Total reads"),
    /** Reads that appear to come from  entirely within a single contig */
    INTERNAL_READS("Single contig reads"),
    /** Number of paired end reads */
    CROSS_CONTIG("Cross contig reads");

    final String mName;
    Stat(String name) {
      mName = name;
    }
  }

  final Map<Stat, Counter> mMap = new EnumMap<>(Stat.class);

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   */
  PacBioStatistics(File outputDirectory) {
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
  public void accumulate(PacBioStatistics stats) {
    for (Map.Entry<Stat, Counter> entry : stats.mMap.entrySet()) {
      mMap.get(entry.getKey()).mCount += entry.getValue().mCount;
    }
  }

  @Override
  public void generateReport() {
  }
}
