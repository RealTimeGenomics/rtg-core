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

import com.rtg.variant.sv.bndeval.BreakpointGeometry;
import com.rtg.variant.sv.bndeval.Orientation;

import junit.framework.TestCase;

/**
 */
public class DebugComparatorTest extends TestCase {
  public void testComparator() {
    final BreakpointConstraint bg1 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 188, 198), 0, 10.0);
    final DiscordantReadSet drs1 = new DiscordantReadSet("f", 60, bg1);

    final DiscordantReadSet drs2 = new DiscordantReadSet("f", 60, bg1);
    assertEquals(0, new DebugComparator().compare(drs1, drs2));

    BreakpointConstraint bg2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "g", 41, 51, 145, 161, 191, 202), 0, 10.0);
    DiscordantReadSet drs3 = new DiscordantReadSet("f", 60, bg2);
    assertEquals(Integer.compare(42, 41), new DebugComparator().compare(drs1, drs3));

    bg2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "g", 42, 51, 145, 148, 187, 199), 0, 10.0);
    drs3 = new DiscordantReadSet("f", 60, bg2);
    assertEquals(Integer.compare(50, 51), new DebugComparator().compare(drs1, drs3));

    bg2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "g", 42, 50, 145, 148, 187, 198), 0, 10.0);
    drs3 = new DiscordantReadSet("f", 60, bg2);
    assertEquals("s".compareTo("g"), new DebugComparator().compare(drs1, drs3));

    bg2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 145, 148, 187, 198), 0, 10.0);
    drs3 = new DiscordantReadSet("f", 60, bg2);
    assertEquals(Integer.compare(142, 145), new DebugComparator().compare(drs1, drs3));

    bg2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 148, 187, 198), 0, 10.0);
    drs3 = new DiscordantReadSet("f", 60, bg2);
    assertEquals(Integer.compare(150, 148), new DebugComparator().compare(drs1, drs3));

    bg2 = new BreakpointConstraint(new BreakpointGeometry(Orientation.UU, "f", "s", 42, 50, 142, 150, 188, 198), 0, 10.0);
    drs3 = new DiscordantReadSet("f", 60, bg2);
    drs3.add(bg2);
    drs3.add(bg2);
    assertEquals(Integer.compare(1, 3), new DebugComparator().compare(drs1, drs3));
  }
}
