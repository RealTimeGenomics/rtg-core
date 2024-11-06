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
import com.rtg.util.License;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;


/**
 */
public class PairedEndMapStatistics extends AbstractStatistics implements MapStatistics {
  private static final String ARM_HEADER = "ARM MAPPINGS";
  private static final String READ_HEADER = "READ MAPPINGS";

  // Column names in header
  private static final String HEADER_LEFT = "left";
  private static final String HEADER_RIGHT = "right";
  private static final String HEADER_BOTH = "both";
  static final int MAX_HEADER_LENGTH = 5; // max length of above columns

  private final SingleEndMapStatistics mLeft;
  private final SingleEndMapStatistics mRight;

  private final boolean mWriteUnmated;

  private int mBothUnmapped = 0;

  /**
   * Constructor
   * @param writeUnmated whether paired-end statistics should include unmated results.
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   */
  public PairedEndMapStatistics(boolean writeUnmated, File outputDirectory) {
    super(outputDirectory);
    mWriteUnmated = writeUnmated;
    mLeft = new SingleEndMapStatistics(outputDirectory);
    mRight = new SingleEndMapStatistics(outputDirectory);
  }

  @Override
  public void increment(MapStatisticsField field, Arm arm) {
    switch (arm) {
    case RIGHT: mRight.increment(field, arm); break;
    case LEFT:
    default:
      mLeft.increment(field, arm);
      break;
    }
  }

  /**
   * Update the special both arms unmapped count.
   */
  public void incrementBothUnmapped() {
    ++mBothUnmapped;
  }

  private static String percentString(double num) {
    return Utils.realFormat(num, 1);
  }

  void appendValue(StringBuilder sb, MapStatisticsField key, String msg, int formatLength) {
    appendValue(sb, key, msg, formatLength, false);
  }

  void appendValue(StringBuilder sb, MapStatisticsField key, String msg, int formatLength, boolean skipIfZero) {
    final long leftValue = value(key, Arm.LEFT);
    final long rightValue = value(key, Arm.RIGHT);
    if (skipIfZero && leftValue == 0 && rightValue == 0) {
      return;
    }
    final long bothValue = totalValue(key);
    final double percent = totalValueAsPercent(key);

    String formatStr = "%" + formatLength + "d";
    final String leftValueFormat = String.format(formatStr, leftValue);
    final String rightValueFormat = String.format(formatStr, rightValue);
    final String bothValueFormat = String.format(formatStr, bothValue);
    formatStr = "%5s";
    final String percentFormat = String.format(formatStr, percentString(percent));
    final String str = leftValueFormat + " " + rightValueFormat + " " + bothValueFormat + " " + percentFormat + "% " + msg;
    sb.append(str).append(StringUtils.LS);
  }

  void appendHeader(StringBuilder sb, int formatLength) {
    assert formatLength >= MAX_HEADER_LENGTH;
    sb.append(ARM_HEADER).append(StringUtils.LS);

    final StringBuilder tsb = new StringBuilder();
    for (int i = HEADER_LEFT.length(); i < formatLength; ++i) {
      tsb.append(" ");
    }
    tsb.append(HEADER_LEFT);
    for (int i = HEADER_RIGHT.length(); i <= formatLength; ++i) {
      tsb.append(" ");
    }
    tsb.append(HEADER_RIGHT);
    for (int i = HEADER_BOTH.length(); i <= formatLength; ++i) {
      tsb.append(" ");
    }
    tsb.append(HEADER_BOTH);
    sb.append(tsb.toString()).append(StringUtils.LS);
  }

  @Override
  protected String getStatistics() {
    return getStatistics(true);
  }

  protected String getStatistics(boolean includeDevLog) {
    final long totalReads = totalValue(MapStatisticsField.TOTAL_READS);
    final int formatLength = Math.max(MAX_HEADER_LENGTH, String.format("%d", totalReads).length());
    final StringBuilder sb = new StringBuilder(StringUtils.LS); //adding extra newLine in case Progress is on
    final StringBuilder devsb = new StringBuilder(StringUtils.LS);

    // header
    appendHeader(sb, formatLength);

    // mated/mapped stats
    appendValue(sb, MapStatisticsField.MATED_UNIQUE_READS, "mated uniquely (NH = 1)", formatLength);
    appendValue(sb, MapStatisticsField.MATED_AMBIG_READS, "mated ambiguously (NH > 1)", formatLength);

    if (mWriteUnmated) {
      appendValue(sb, MapStatisticsField.UNMATED_UNIQUE_READS, "unmated uniquely (NH = 1)", formatLength);
      appendValue(sb, MapStatisticsField.UNMATED_AMBIG_READS, "unmated ambiguously (NH > 1)", formatLength);
    }

    //unmapped stats
    appendValue(sb, MapStatisticsField.UNMAPPED_BLOCKED, "unmapped due to read frequency (XC = B)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_TOPN, "unmapped with no matings but too many hits (XC = C)" , formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_MATED_POOR, "unmapped with poor matings (XC = d)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_MATED_TOO_MANY, "unmapped with too many matings (XC = e)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_UNMATED_POOR, "unmapped with no matings and poor hits (XC = D)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_UNMATED_TOO_MANY, "unmapped with no matings and too many good hits (XC = E)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_NO_HITS, "unmapped with no hits (XC = A)", formatLength);

    // always print total
    appendValue(sb, MapStatisticsField.TOTAL_READS, "total", formatLength);

    if (License.isDeveloper()) {
      appendValue(sb, MapStatisticsField.MISSING, "arms missing", formatLength, true);
      appendValue(devsb, MapStatisticsField.IGNORED, "arms ignored", formatLength, true);

      final long matedReads = totalValue(MapStatisticsField.MATED_AMBIG_READS)
      + value(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT);
      final long mappedReads = totalReads - mBothUnmapped;
      final long unmatedReads = mappedReads - matedReads;
      final int formatLength2 = String.format("%d", totalReads).length();
      devsb.append(READ_HEADER).append(StringUtils.LS);
      appendValue(devsb, totalReads, totalReads, "total reads", formatLength2);
      appendValue(devsb, matedReads, totalReads, "mated reads", formatLength2);
      appendValue(devsb, unmatedReads, totalReads, "unmated reads", formatLength2);
      appendValue(devsb, mappedReads, totalReads, "mapped reads", formatLength2);
      appendValue(devsb, mBothUnmapped, totalReads, "unmapped reads", formatLength2);

      Diagnostic.developerLog(devsb.toString());
    }
    return sb.toString();
  }

  private static String percent(long numerator, long denominator) {
    return percentString(100.0 * numerator / denominator);
  }

  private void appendValue(StringBuilder sb, long value, long total, String msg, int formatLength) {
    String formatStr = "%" + formatLength + "d";
    final String valueFormat = String.format(formatStr, value);
    formatStr = "%5s";
    final String percentFormat = String.format(formatStr, percent(value, total));
    final String str = valueFormat + " " + percentFormat + "% " + msg;
    sb.append(str).append(StringUtils.LS);
  }

  @Override
  public void reset() {
    mLeft.reset();
    mRight.reset();
    mBothUnmapped = 0;
  }

  @Override
  public long value(MapStatisticsField field, Arm arm) {
    final long res;
    switch (arm) {
      case RIGHT:
        res = mRight.value(field, arm);
        break;
      case LEFT:
      default:
        res = mLeft.value(field, arm);
        break;
    }
    return res;
  }

  @Override
  public long totalValue(MapStatisticsField field) {
    return mLeft.totalValue(field) + mRight.totalValue(field);
  }

  @Override
  public double valueAsPercent(MapStatisticsField field, Arm arm) {
    final double res;
    switch (arm) {
    case RIGHT: res = mRight.valueAsPercent(field, arm); break;
    case LEFT:
    default:
      res = mLeft.valueAsPercent(field, arm); break;
    }
    return res;
  }

  @Override
  public double totalValueAsPercent(MapStatisticsField field) {
    final long total = mLeft.totalValue(MapStatisticsField.TOTAL_READS) + mRight.totalValue(MapStatisticsField.TOTAL_READS);
     return (mLeft.totalValue(field) + mRight.totalValue(field)) * 100.0 / (total <= 0 ? 1 : total);
  }

  @Override
  public void set(MapStatisticsField field, Arm arm, long value) {
    switch (arm) {
    case RIGHT: mRight.set(field, arm, value); break;
    case LEFT:
    default:
      mLeft.set(field, arm, value);
      break;
    }
  }

  @Override
  public void generateReport() {
  }
}
