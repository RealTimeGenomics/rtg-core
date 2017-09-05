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

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import com.rtg.util.Environment;
import com.rtg.util.cli.CommandLine;

/**
 * Produces simple BED output for breakends
 */
public class BedDiscordantOutputFormatter {

  /**
   * Format the given object
   * @param readset read set to output
   * @return string containing formatted results
   */
  public DiscordBedRecord format(DiscordantReadSet readset) {
    final BreakpointConstraint geo = readset.getIntersection() == null ? readset.getUnion() : readset.getIntersection();
    final BreakpointPosition pos = geo.position(); // Do magic
    final BreakpointPosition pos2 = geo.flip().position(); // Do magic on the other end

    //Use Math.max(x, 0) to truncate records that go before the start of the reference
    final DiscordBedRecord record = new DiscordBedRecord(readset.getSequenceName(), Math.max(pos.lo(), 0), Math.max(pos.hi(), 0),
        "remote:" + geo.getYName() + ":" + Math.max(pos2.lo(), 0) + "-" + Math.max(pos2.hi(), 0), // Name
        "" + readset.getCounts()); // Score
    if (readset.getIntersection() == null) {
      record.setFiltered();
    }
    return record;
  }

  /**
   * @return header for current format
   */
  public String header() {
    final StringBuilder sb = new StringBuilder();
    sb.append("#Version ").append(Environment.getVersion()).append(", Discordance output ").append(DiscordantTool.DT_OUTPUT_VERSION).append(LS);
    if (CommandLine.getCommandLine() != null) {
      sb.append("#CL" + TAB).append(CommandLine.getCommandLine()).append(LS);
    }
    sb.append("#RUN-ID" + TAB).append(CommandLine.getRunId()).append(LS);
    sb.append("#");
    sb.append("chromosome" + TAB + "start" + TAB + "end" + TAB);
    sb.append("remote" + TAB + "count");
    sb.append(LS);
    return sb.toString();
  }

}
