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

import com.rtg.reader.Arm;
import com.rtg.util.StringUtils;

/**
 * Statistics for map filter
 */
public class MapFilterPairedMapStatistics extends PairedEndMapStatistics {

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   */
  public MapFilterPairedMapStatistics(File outputDirectory) {
    super(true, outputDirectory);
  }

  private void checkField(MapStatisticsField field) {
    switch (field) {
      case UNMATED_AMBIG_READS:
      case UNMAPPED_UNMATED_POOR:
      case UNMAPPED_NO_HITS:
      case TOTAL_READS:
        break;
      default:
        throw new UnsupportedOperationException("Field " + field + " is not supported for this command");
    }
  }

  @Override
  public void increment(MapStatisticsField field, Arm arm) {
    checkField(field);
    super.increment(field, arm);
  }

  @Override
  public void set(MapStatisticsField field, Arm arm, long value) {
    checkField(field);
    super.set(field, arm, value);
  }

  @Override
  public long value(MapStatisticsField field, Arm arm) {
    checkField(field);
    return super.value(field, arm);
  }

  @Override
  public long totalValue(MapStatisticsField field) {
    checkField(field);
    return super.totalValue(field);
  }

  @Override
  public double valueAsPercent(MapStatisticsField field, Arm arm) {
    checkField(field);
    return super.valueAsPercent(field, arm);
  }

  @Override
  public double totalValueAsPercent(MapStatisticsField field) {
    checkField(field);
    return super.totalValueAsPercent(field);
  }


  @Override
  protected String getStatistics(boolean includeDevLog) {
    final long total = totalValue(MapStatisticsField.TOTAL_READS);
    final int formatLength = Math.max(MAX_HEADER_LENGTH, String.format("%d", total).length());
    final StringBuilder sb = new StringBuilder();
    //adding extra newLine in case Progress is on
    sb.append(StringUtils.LS);

    // header
    appendHeader(sb, formatLength);

    // mated/mapped stats
    appendValue(sb, MapStatisticsField.UNMATED_AMBIG_READS, "mapped", formatLength);
    //unmapped stats
    appendValue(sb, MapStatisticsField.UNMAPPED_UNMATED_POOR, "unmapped with poor hits (XC = D)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_NO_HITS, "unmapped with no hits (XC = A)", formatLength);
    // always print total
    appendValue(sb, MapStatisticsField.TOTAL_READS, "total", formatLength);
    return sb.toString();
  }
}
