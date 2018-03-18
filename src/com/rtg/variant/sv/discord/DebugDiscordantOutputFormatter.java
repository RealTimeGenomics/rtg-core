/*
 * Copyright (c) 2018. Real Time Genomics Limited.
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

import java.io.IOException;
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
    public BlockCompressedPositionReader getReader(InputStream is) throws IOException {
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
