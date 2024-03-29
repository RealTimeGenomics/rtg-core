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
package com.rtg.ngs;

import java.io.File;

import com.rtg.reader.Arm;
import com.rtg.util.StringUtils;

/**
 * Statistics for map filter
 */
public class MapFilterSingleEndMapStatistics extends SingleEndMapStatistics {

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   */
  public MapFilterSingleEndMapStatistics(File outputDirectory) {
    super(outputDirectory);
  }

  private void checkField(MapStatisticsField field) {
     switch (field) {
      case UNMATED_AMBIG_READS:
      case UNMAPPED_NO_HITS:
      case UNMAPPED_UNMATED_POOR:
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
  protected String getStatistics() {
    final int formatLength = String.format("%d", mTotal).length();
    final StringBuilder sb = new StringBuilder();
    //adding extra newLine in case Progress is on
    sb.append(StringUtils.LS);

    // header lines
    sb.append(HEADER).append(StringUtils.LS);

    // mapped stats
    appendValue(sb, MapStatisticsField.UNMATED_AMBIG_READS, "mapped", formatLength);

    //unmapped stats
    appendValue(sb, MapStatisticsField.UNMAPPED_UNMATED_POOR, "unmapped with poor hits (XC = D)", formatLength);
    appendValue(sb, MapStatisticsField.UNMAPPED_NO_HITS, "unmapped with no hits (XC = A)", formatLength);

    appendValue(sb, MapStatisticsField.TOTAL_READS, "total", formatLength);
    return sb.toString();
  }

}
