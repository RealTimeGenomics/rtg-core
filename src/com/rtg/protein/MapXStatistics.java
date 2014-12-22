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

  @Override
  public void generateReport() {
    //TODO somehow wrangle the MapXCli report stuff into here - unfortunately it heavily relies on the params object of all things, so that'll be fun.
  }
}
