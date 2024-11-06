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
