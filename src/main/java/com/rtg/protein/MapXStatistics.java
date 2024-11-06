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
package com.rtg.protein;

import java.io.File;

import com.rtg.ngs.MapStatisticsField;
import com.rtg.ngs.SingleEndMapStatistics;
import com.rtg.util.StringUtils;


/**
 */
public class MapXStatistics extends SingleEndMapStatistics {

  private static final String MAPX_HEADER = "ALIGNMENT STATISTICS";

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   */
  public MapXStatistics(File outputDirectory) {
    super(outputDirectory);
  }

  @Override
  protected String getStatistics() {
    final int formatLength = String.format("%d", mTotal).length();
    final StringBuilder sb = new StringBuilder();
    //adding extra newLine in case Progress is on
    sb.append(StringUtils.LS);

    // header lines
    sb.append(MAPX_HEADER).append(StringUtils.LS);

    appendValue(sb, MapStatisticsField.UNMATED_UNIQUE_READS, "alignments met threshold", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_UNMATED_POOR, "alignments failed threshold", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_NO_HITS, "no significant alignments", formatLength);

    appendValue(sb, MapStatisticsField.TOTAL_READS, "total", formatLength);
    return sb.toString();
  }
}
