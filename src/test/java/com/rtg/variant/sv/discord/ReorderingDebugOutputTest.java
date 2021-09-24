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
