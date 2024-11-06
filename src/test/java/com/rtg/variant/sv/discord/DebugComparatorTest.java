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
