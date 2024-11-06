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

import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.sv.bndeval.BreakpointGeometry;
import com.rtg.variant.sv.bndeval.Orientation;

import junit.framework.TestCase;

/**
 */
public class BedDiscordantOutputFormatterTest extends TestCase {

  public BedDiscordantOutputFormatterTest(String name) {
    super(name);
  }

  public void test() {
    final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
    final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);
    assertEquals(42, drs.unionPosition());
    assertEquals(102, drs.flushPosition());
    final DiscordBedRecord rec = new BedDiscordantOutputFormatter().format(drs);

    assertEquals("f", rec.getSequenceName());
    assertEquals(42, rec.getStart());
    assertEquals(48, rec.getEnd());
    assertEquals(2, rec.getAnnotations().length);
    assertFalse(rec.isFiltered());
  }

  public void testNegative() throws IOException {
    Diagnostic.setLogStream();
    final File t = FileUtils.createTempDir("vcfdescord", "test");
    try {
      final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "f", -20, 2, -40, 3, -60, -17), 0, 10.0); //new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
      final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);
      assertEquals(-20, drs.unionPosition());
      assertEquals(40, drs.flushPosition());
      final String exp = "f" + TAB
          + "0" + TAB
          + "2" + TAB
          + "remote:f:0-2" + TAB
          + "1";
      assertEquals(exp, new BedDiscordantOutputFormatter().format(drs).toString());
    } finally {
      assertTrue(FileHelper.deleteAll(t));
    }
  }

  public void testNegative2() throws IOException {
    Diagnostic.setLogStream();
    final File t = FileUtils.createTempDir("vcfdescord", "test");
    try {
      final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.DD, "f", "f", 2, -20, 3, -40, 0, 55), 0, 10.0); //new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
      final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);
      assertEquals(2, drs.unionPosition());
      assertEquals(62, drs.flushPosition());
      final String exp = "f" + TAB
          + "0" + TAB
          + "2" + TAB
          + "remote:f:0-3" + TAB
          + "1";
      assertEquals(exp, new BedDiscordantOutputFormatter().format(drs).toString());
    } finally {
      assertTrue(FileHelper.deleteAll(t));
    }
  }

  public void testFilter() {
    final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
    final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);
    final BreakpointConstraint bg2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 842, 850, 942, 950, 1792, 1796), 550, 50.0);
    final DiscordantReadSet drs2 = new DiscordantReadSet("f", 860, bg2);
    drs.addAll(drs2);

    final DiscordBedRecord rec = new BedDiscordantOutputFormatter().format(drs);
    assertTrue(rec.isFiltered());
  }
}
