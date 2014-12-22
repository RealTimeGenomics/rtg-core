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

import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;


/**
 */
public class DiscordantToolStatisticsTest extends TestCase {

  private NanoRegression mNano = null;

  @Override
  public void setUp() {
    mNano = new NanoRegression(this.getClass());
  }

  @Override
  public void tearDown() throws IOException {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  public void testStatisticsFiltered() throws IOException {
    final DiscordantToolStatistics stats = new DiscordantToolStatistics(null);

    stats.tallyDiscordantReadSet(new DiscordantReadSet("Blah", 1, new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 1.0, 0.5)));
    stats.tallyDiscordantReadSet(new DiscordantReadSet("Blah", 1, new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 1.0, 0.5)));

    final BreakpointConstraint bg = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 190, 198), 0, 10.0);
    final DiscordantReadSet drs = new DiscordantReadSet("f", 60, bg);
    final BreakpointConstraint bg2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 842, 850, 942, 950, 1792, 1796), 550, 50.0);
    final DiscordantReadSet drs2 = new DiscordantReadSet("f", 860, bg2);
    drs.addAll(drs2);
    stats.tallyDiscordantReadSet(drs);
    mNano.check("discordantstats", stats.getStatistics());
  }
}
