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

import static com.rtg.util.StringUtils.TAB;

import java.io.File;
import java.io.IOException;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

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
