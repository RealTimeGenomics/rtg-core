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

import com.rtg.util.TestUtils;
import com.rtg.variant.sv.bndeval.BreakpointGeometry;
import com.rtg.variant.sv.bndeval.Orientation;

import junit.framework.TestCase;

/**
 */
public class DebugDiscordantOutputFormatterTest extends TestCase {

  public DebugDiscordantOutputFormatterTest(String name) {
    super(name);
  }

  public void testHeader() {
    final String str = new DebugDiscordantOutputFormatter().header().replace("\t", " ");
    TestUtils.containsAll(str, "#Version ",
      "Discordance output 1",
      "#RUN-ID ",
      "#count union:template start end remote start end r s intersection:template start end remote start end r s");
  }

  public void test() {
    final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
    final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);
    assertEquals(42, drs.unionPosition());
    assertEquals(102, drs.flushPosition());
    assertEquals(""
        + "1" + TAB + "f" + TAB + "42" + TAB + "50" + TAB + "s" + TAB + "142" + TAB + "150" + TAB + "190" + TAB + "198" + TAB
        + "f" + TAB + "42" + TAB + "50" + TAB + "s" + TAB + "142" + TAB + "150" + TAB + "190" + TAB + "198" + LS
        , new DebugDiscordantOutputFormatter().format(drs));
  }

}
