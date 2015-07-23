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
 * Produces debug output for investigations
 */
public class DebugDiscordantOutputFormatter {

  /**
   * Format the given object
   * @param readset read set to output
   * @return string containing formatted results
   */
  public String format(DiscordantReadSet readset) {
    final StringBuilder sb = new StringBuilder();
    sb.append(readset.getCounts());
    sb.append(TAB);
    output(readset.getUnion(), sb);
    if (readset.getIntersection() != null) {
      sb.append(TAB);
      output(readset.getIntersection(), sb);
    }
    sb.append(LS);
    return sb.toString();
  }

  private void output(AbstractBreakpointGeometry geo, StringBuilder sb) {
    sb.append(geo.getXName());
    sb.append(TAB);
    sb.append(geo.getX());
    sb.append(TAB);
    sb.append(geo.getZ());
    sb.append(TAB);
    sb.append(geo.getYName());
    sb.append(TAB);
    sb.append(geo.getY());
    sb.append(TAB);
    sb.append(geo.getW());
    sb.append(TAB);
    sb.append(geo.getR());
    sb.append(TAB);
    sb.append(geo.getS());
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
    sb.append("count" + TAB);
    sb.append("union:template" + TAB + "start" + TAB + "end" + TAB);
    sb.append("remote" + TAB + "start" + TAB + "end" + TAB);
    sb.append("r" + TAB + "s" + TAB);
    sb.append("intersection:template" + TAB + "start" + TAB + "end" + TAB);
    sb.append("remote" + TAB + "start" + TAB + "end" + TAB);
    sb.append("r" + TAB + "s" + TAB);
    sb.append(LS);
    return sb.toString();
  }

}
