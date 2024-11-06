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
