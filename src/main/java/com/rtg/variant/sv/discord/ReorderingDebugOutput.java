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
