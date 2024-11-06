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
  public void accumulate(PacBioStatistics stats) {
    for (Map.Entry<Stat, Counter> entry : stats.mMap.entrySet()) {
      mMap.get(entry.getKey()).mCount += entry.getValue().mCount;
    }
  }

  @Override
  public void generateReport() {
  }
}
