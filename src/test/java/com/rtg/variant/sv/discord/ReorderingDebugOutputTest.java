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

import java.io.IOException;

import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.variant.sv.bndeval.BreakpointGeometry;
import com.rtg.variant.sv.bndeval.Orientation;

import junit.framework.TestCase;

/**
 */
public class ReorderingDebugOutputTest extends TestCase {
  public void test() throws IOException {
    BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
    final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);

    bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 43, 50, 142, 150, 190, 198), 0, 10.0);
    final DiscordantReadSet drs2 = new DiscordantReadSet("f", 60, bg);

    final MemoryPrintStream mps = new MemoryPrintStream();
    try (ReorderingDebugOutput debug = new ReorderingDebugOutput(new DebugDiscordantOutputFormatter(), mps.outputStream(), 100)) {
      debug.addRecord(drs2);
      debug.addRecord(drs);
    }
    final String expected = ("1 f 42 50 s 142 150 190 198 f 42 50 s 142 150 190 198" + LS
     + "1 f 43 50 s 142 150 190 198 f 43 50 s 142 150 190 198" + LS
    ).replace(' ', '\t');
    assertEquals(expected, mps.toString());

  }

  public void testLimit() throws IOException {
    final MemoryPrintStream diag = new MemoryPrintStream();
    Diagnostic.setLogStream(diag.printStream());
    try {
      BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 184, 192), 0, 10.0);
      final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);

      bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 143, 150, 142, 150, 285, 293), 0, 10.0);
      final DiscordantReadSet drs2 = new DiscordantReadSet("f", 60, bg);
      bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 41, 150, 142, 150, 190, 292), 0, 10.0);
      final DiscordantReadSet drs3 = new DiscordantReadSet("f", 60, bg);

      final MemoryPrintStream mps = new MemoryPrintStream();
      try (ReorderingDebugOutput debug = new ReorderingDebugOutput(new DebugDiscordantOutputFormatter(), mps.outputStream(), 100)) {
        debug.addRecord(drs2);
        debug.addRecord(drs);
        debug.addRecord(drs3);

      }
    } finally {
      Diagnostic.setLogStream();
    }
    TestUtils.containsAll(diag.toString(), "Debug breakpoint output dropped due to reordering failure");
  }
}
