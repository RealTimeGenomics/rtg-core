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
        String.valueOf(readset.getCounts())); // Score
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
