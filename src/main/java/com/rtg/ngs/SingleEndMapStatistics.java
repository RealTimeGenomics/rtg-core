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

import java.io.File;

import com.rtg.launcher.AbstractStatistics;
import com.rtg.reader.Arm;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;

/**
 */
public class SingleEndMapStatistics extends AbstractStatistics implements MapStatistics {
  protected static final String HEADER = "READ MAPPINGS";

  private long mMissing = 0;
  private long mIgnored = 0;
  private long mMatedUnique = 0;
  private long mMatedAmbig = 0;
  private long mUnmatedUnique = 0;
  private long mUnmatedAmbig = 0;
  private long mUnmappedBlocked = 0;      //XC B
  private long mUnmappedMatedTooMany = 0; //XC e
  private long mUnmappedMatedPoor = 0;    //XC d
  private long mUnmappedUnmatedTooMany = 0;   //XC E
  private long mUnmappedUnmatedPoor = 0;  //XC D
  private long mUnmappedTopN = 0;  //XC C
  private long mUnmappedNoHits = 0;
  protected long mTotal = 0;

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   */
  public SingleEndMapStatistics(File outputDirectory) {
    super(outputDirectory);
  }

  @Override
  public void increment(MapStatisticsField field, Arm arm) {
    switch (field) {
    case MATED_UNIQUE_READS: mMatedUnique++; break;
    case MATED_AMBIG_READS: mMatedAmbig++; break;
    case UNMATED_UNIQUE_READS: mUnmatedUnique++; break;
    case UNMATED_AMBIG_READS: mUnmatedAmbig++; break;
    case UNMAPPED_NO_HITS: mUnmappedNoHits++; break;
    case UNMAPPED_BLOCKED: mUnmappedBlocked++; break;
    case UNMAPPED_MATED_POOR: mUnmappedMatedPoor++; break;
    case UNMAPPED_MATED_TOO_MANY: mUnmappedMatedTooMany++; break;
    case UNMAPPED_TOPN: mUnmappedTopN++; break;
    case UNMAPPED_UNMATED_POOR: mUnmappedUnmatedPoor++; break;
    case UNMAPPED_UNMATED_TOO_MANY: mUnmappedUnmatedTooMany++; break;
    case MISSING: mMissing++ ; break;
    case IGNORED: mIgnored++; break;
    case TOTAL_READS: mTotal++; break;
    default:
      throw new RuntimeException();
    }
  }

  private static String percentString(double num) {
    return Utils.realFormat(num, 1);
  }

  protected void appendValue(StringBuilder sb, MapStatisticsField key, String msg, int formatLength) {
    final long value = totalValue(key);
    final double percent = totalValueAsPercent(key);

    String formatStr = "%" + formatLength + "d";
    final String valueFormat = String.format(formatStr, value);
    formatStr = "%5s";
    final String percentFormat = String.format(formatStr, percentString(percent));
    final String str = valueFormat + " " + percentFormat + "% " + msg;
    sb.append(str).append(StringUtils.LS);
  }

  @Override
  protected String getStatistics() {
    final int formatLength = String.format("%d", mTotal).length();
    final StringBuilder sb = new StringBuilder();
    //adding extra newLine in case Progress is on
    sb.append(StringUtils.LS);

    // header lines
    sb.append(HEADER).append(StringUtils.LS);

    // mapped stats
    appendValue(sb, MapStatisticsField.UNMATED_UNIQUE_READS, "mapped uniquely (NH = 1)", formatLength);
    appendValue(sb, MapStatisticsField.UNMATED_AMBIG_READS, "mapped ambiguously (NH > 1)", formatLength);

    //unmapped stats
    appendValue(sb, MapStatisticsField.UNMAPPED_BLOCKED, "unmapped due to read frequency (XC = B)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_TOPN, "unmapped with too many hits (XC = C)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_UNMATED_POOR, "unmapped with poor hits (XC = D)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_UNMATED_TOO_MANY, "unmapped with too many good hits (XC = E)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_NO_HITS, "unmapped with no hits (XC = A)", formatLength);

    appendValue(sb, MapStatisticsField.TOTAL_READS, "total", formatLength);
    return sb.toString();
  }

  @Override
  public void reset() {
    mUnmappedNoHits = 0;
    mMissing = 0;
    mIgnored = 0;
    mMatedUnique = 0;
    mMatedAmbig = 0;
    mUnmatedUnique = 0;
    mUnmatedAmbig = 0;
    mUnmappedBlocked = 0;      //XC B
    mUnmappedMatedTooMany = 0; //XC e
    mUnmappedMatedPoor = 0;    //XC d
    mUnmappedUnmatedTooMany = 0;   //XC E
    mUnmappedUnmatedPoor = 0;  //XC D
    mUnmappedTopN = 0;  //XC C
    mTotal = 0;
  }

  @Override
  public long value(MapStatisticsField field, Arm arm) {
    return totalValue(field);
  }

  @Override
  public long totalValue(MapStatisticsField field) {
    final long res;
    switch (field) {
    case MATED_UNIQUE_READS: res = mMatedUnique; break;
    case MATED_AMBIG_READS: res = mMatedAmbig; break;
    case UNMATED_UNIQUE_READS: res = mUnmatedUnique; break;
    case UNMATED_AMBIG_READS: res = mUnmatedAmbig; break;
    case UNMAPPED_NO_HITS: res = mUnmappedNoHits; break;
    case UNMAPPED_BLOCKED: res = mUnmappedBlocked; break;
    case UNMAPPED_MATED_POOR: res = mUnmappedMatedPoor; break;
    case UNMAPPED_MATED_TOO_MANY: res = mUnmappedMatedTooMany; break;
    case UNMAPPED_TOPN: res = mUnmappedTopN; break;
    case UNMAPPED_UNMATED_POOR: res = mUnmappedUnmatedPoor; break;
    case UNMAPPED_UNMATED_TOO_MANY: res = mUnmappedUnmatedTooMany; break;
    case MISSING: res = mMissing; break;
    case IGNORED: res = mIgnored; break;
    case TOTAL_READS: res = mTotal; break;
    default:
      res = -1;
    }
    return res;
  }

  @Override
  public double valueAsPercent(MapStatisticsField field, Arm arm) {
    return totalValueAsPercent(field);
  }

  @Override
  public double totalValueAsPercent(MapStatisticsField field) {
    return totalValue(field) * 100.0 / ((mTotal <= 0) ? 1 : mTotal);
  }

  @Override
  public void set(MapStatisticsField field, Arm arm, long value) {
    switch (field) {
    case MATED_UNIQUE_READS: mMatedUnique = value; break;
    case MATED_AMBIG_READS: mMatedAmbig = value; break;
    case UNMATED_UNIQUE_READS: mUnmatedUnique = value; break;
    case UNMATED_AMBIG_READS: mUnmatedAmbig = value; break;
    case UNMAPPED_NO_HITS: mUnmappedNoHits = value; break;
    case UNMAPPED_BLOCKED: mUnmappedBlocked = value; break;
    case UNMAPPED_MATED_POOR: mUnmappedMatedPoor = value; break;
    case UNMAPPED_MATED_TOO_MANY: mUnmappedMatedTooMany = value; break;
    case UNMAPPED_TOPN: mUnmappedTopN = value; break;
    case UNMAPPED_UNMATED_POOR: mUnmappedUnmatedPoor = value; break;
    case UNMAPPED_UNMATED_TOO_MANY: mUnmappedUnmatedTooMany = value; break;
    case MISSING: mMissing = value ; break;
    case TOTAL_READS: mTotal = value; break;
    default:
      throw new InternalError();
    }
  }

  @Override
  public void generateReport() {
  }
}
