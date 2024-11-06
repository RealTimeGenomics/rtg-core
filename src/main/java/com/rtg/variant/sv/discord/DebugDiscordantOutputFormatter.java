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

import java.io.InputStream;

import com.rtg.tabix.BlockCompressedLineReader;
import com.rtg.tabix.BlockCompressedPositionReader;
import com.rtg.tabix.GenericPositionReader;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.Environment;
import com.rtg.util.cli.CommandLine;
import com.rtg.variant.sv.bndeval.AbstractBreakpointGeometry;

import htsjdk.samtools.util.BlockCompressedInputStream;

/**
 * Produces debug output for investigations
 */
public class DebugDiscordantOutputFormatter {

  /** Permit tabix indexing of debug output files */
  static class DebugIndexerFactory extends TabixIndexer.IndexerFactory {

    DebugIndexerFactory() {
      super(0);
    }

    @Override
    public TabixIndexer.TabixOptions getOptions() {
      return new TabixIndexer.TabixOptions(TabixIndexer.TabixOptions.FORMAT_GENERIC, 1, 2, 3, '#', mSkip, true);
    }

    @Override
    public BlockCompressedPositionReader getReader(InputStream is) {
      return new GenericPositionReader(new BlockCompressedLineReader(new BlockCompressedInputStream(is)), getOptions());
    }
  }

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
    sb.append(geo.getXLo());
    sb.append(TAB);
    sb.append(geo.getXHi());
    sb.append(TAB);
    sb.append(geo.getYName());
    sb.append(TAB);
    sb.append(geo.getYLo());
    sb.append(TAB);
    sb.append(geo.getYHi());
    sb.append(TAB);
    sb.append(geo.getRLo());
    sb.append(TAB);
    sb.append(geo.getRHi());
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
