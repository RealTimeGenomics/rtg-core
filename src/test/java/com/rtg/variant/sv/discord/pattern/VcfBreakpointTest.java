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

package com.rtg.variant.sv.discord.pattern;

import com.rtg.vcf.VcfRecord;

import junit.framework.TestCase;

public class VcfBreakpointTest extends TestCase {
  public void test() {
    final VcfBreakpoint breakpoint = makeBreakpoint("asdf", 101, "A[chr1:301[");
    assertEquals("VcfBreakpoint: asdf 100 A[chr1:301[", breakpoint.toString());
    assertEquals("asdf", breakpoint.getLocalChr());
    assertEquals(100, breakpoint.getLocalPos());
    assertEquals("chr1", breakpoint.getRemoteChr());
    assertEquals(300, breakpoint.getRemotePos());
    assertEquals(0, breakpoint.getDepth());
    assertTrue(breakpoint.isLocalUp());
    assertFalse(breakpoint.isRemoteUp());
    VcfBreakpoint breakpoint2 = makeBreakpoint("asdf", 101, "A[chr1:300[");
    assertEquals(0, breakpoint.compareTo(breakpoint));
    assertEquals(1, breakpoint.compareTo(breakpoint2));
    assertEquals(-1, breakpoint2.compareTo(breakpoint));
    breakpoint2 = makeBreakpoint("asde", 101, "A[chr1:301[");
    assertEquals(1, breakpoint.compareTo(breakpoint2));
    breakpoint2 = makeBreakpoint("asdf", 100, "A[chr1:301[");
    assertEquals(1, breakpoint.compareTo(breakpoint2));
    breakpoint2 = makeBreakpoint("asdf", 101, "A[chm:301[");
    assertEquals("chr1".compareTo("chm"), breakpoint.compareTo(breakpoint2));
    breakpoint2 = makeBreakpoint("asdf", 101, "[chr1:301[A");
    assertEquals(1, breakpoint.compareTo(breakpoint2));
    breakpoint2 = makeBreakpoint("asdf", 101, "A]chr1:301]");
    assertEquals(-1, breakpoint.compareTo(breakpoint2));
  }

  private static VcfBreakpoint makeBreakpoint(String local, int pos, String alt) {
    final VcfRecord rec = makeRecord(local, pos, alt);
    return new VcfBreakpoint(rec);
  }

  private static VcfRecord makeRecord(String local, int pos, String alt) {
    final VcfRecord bar = new VcfRecord(local, pos - 1, "A");
    bar.addAltCall(alt);
    return bar;
  }

  public void testDepth() {
    final VcfRecord rec = makeRecord("foo", 20, "A[chr:200[");
    rec.setInfo("DP", "7");
    final VcfBreakpoint breakpoint = new VcfBreakpoint(rec);
    assertEquals(7, breakpoint.getDepth());

  }

  static final VcfBreakpoint BREAK_POINT = makeBreakpoint("foo", 21, "A[chr1:201[");
  public void testEquals() {
    final VcfBreakpoint b2 = makeBreakpoint("foo", 22, "A[chr1:201[");
    assertFalse(BREAK_POINT.equals(null));
    assertFalse(BREAK_POINT.equals("monkey"));
    assertTrue(BREAK_POINT.equals(BREAK_POINT));
    assertFalse(BREAK_POINT.equals(b2));
  }

  public void testHashCode() {
    assertEquals(-734628597, BREAK_POINT.hashCode());
  }
}
