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

import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;

public class BreakpointStoreTest extends TestCase {
  public void testAdding() {
    final BreakpointStore bs = new BreakpointStore();
    bs.add(new VcfBreakpoint("bar", 100, "foo", 20, false, false, 0));
    bs.add(new VcfBreakpoint("foo", 20, "bar", 100, false, false, 0));
    bs.add(new VcfBreakpoint("foo", 22, "bar", 100, false, false, 0));
    bs.add(new VcfBreakpoint("foo", 22, "foo", 101, false, false, 0));
    bs.add(new VcfBreakpoint("foo", 22, "foo", 100, false, false, 0));
    final Iterator<VcfBreakpoint> it = bs.iterator();
    VcfBreakpoint b = it.next();
    assertEquals("VcfBreakpoint: bar 100 [foo:21[N", b.toString());
    b = it.next();
    assertEquals("VcfBreakpoint: foo 20 [bar:101[N", b.toString());
    b = it.next();
    assertEquals("VcfBreakpoint: foo 22 [bar:101[N", b.toString());
    b = it.next();
    assertEquals("VcfBreakpoint: foo 22 [foo:101[N", b.toString());
    b = it.next();
    assertEquals("VcfBreakpoint: foo 22 [foo:102[N", b.toString());
    assertFalse(it.hasNext());
    b = bs.getMap().get("foo").get("bar").first();
    assertEquals("VcfBreakpoint: foo 20 [bar:101[N", b.toString());
    final List<String> names = bs.getChromosomes();
    assertEquals(2, names.size());
    assertEquals("bar", names.get(0));
    assertEquals("foo", names.get(1));
  }
}
