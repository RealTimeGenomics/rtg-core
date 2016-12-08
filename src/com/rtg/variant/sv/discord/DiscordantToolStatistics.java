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
package com.rtg.variant.sv.discord;

import java.io.File;

import com.rtg.launcher.AbstractStatistics;
import com.rtg.util.TextTable;

/**
 * Statistics object for tracking and printing discord modules statistics.
 */
public class DiscordantToolStatistics extends AbstractStatistics {

  protected long mTotalBreakEnds;
  protected long mFiltered;

  /**
   * @param outputDirectory The base output directory to generate statistics and reports in. May be null if no statistics or reports are to be generated.
   */
  public DiscordantToolStatistics(File outputDirectory) {
    super(outputDirectory);
  }

  // For when eventually we do record more statistics about a record set
  protected void tallyDiscordantReadSet(DiscordantReadSet drs) {
    if (drs != null) {
      ++mTotalBreakEnds;
      if (drs.getIntersection() == null) {
        ++mFiltered;
      }
    }
  }

  @Override
  protected String getStatistics() {
    final TextTable table = new TextTable();
    table.addRow("Total calls:", "" + mTotalBreakEnds);
    table.addRow("Passed Filters:", "" + (mTotalBreakEnds - mFiltered));
    table.addRow("Failed Filters:", "" + mFiltered);
    table.addRow("Structural variant breakends:", "" + (mTotalBreakEnds - mFiltered));
    return table.toString();
  }

  @Override
  public void generateReport() {
  }
}
