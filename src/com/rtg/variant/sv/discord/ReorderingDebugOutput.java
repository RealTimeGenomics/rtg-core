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

import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.ReorderingQueue;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Reorders the debug output of the discordant tool since it may be generated in an order that doesn't match desired
 * output order. Duplicate records will be removed
 */
class ReorderingDebugOutput extends ReorderingQueue<DiscordantReadSet> {

  private final OutputStream mDebug;
  private final DebugDiscordantOutputFormatter mFormatter;

  /**
   *
   * @param formatter formats the output
   * @param output the  debug output stream
   * @param bufferDistance the initial buffer distance. If records are out of order by more than this errors will be
   *                       generated and it will be skipped.
   */
  ReorderingDebugOutput(DebugDiscordantOutputFormatter formatter, OutputStream output, int bufferDistance) {
    super(bufferDistance, new DebugComparator());
    mDebug = output;
    mFormatter = formatter;
  }

  @Override
  protected String getReferenceName(DiscordantReadSet record) {
    return record.getSequenceName();
  }

  @Override
  protected int getPosition(DiscordantReadSet record) {
    return record.getUnion().getXLo();
  }

  @Override
  protected void flushRecord(DiscordantReadSet rec) throws IOException {
    final String debugRep = mFormatter.format(rec);
    mDebug.write(debugRep.getBytes());
  }

  @Override
  protected void reportReorderingFailure(DiscordantReadSet rec) {
    Diagnostic.warning("Debug breakpoint output dropped due to reordering failure.");
  }
}
