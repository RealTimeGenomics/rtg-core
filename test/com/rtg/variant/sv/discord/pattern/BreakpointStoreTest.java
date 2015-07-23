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

package com.rtg.variant.sv.discord.pattern;

import java.util.Iterator;
import java.util.List;

import junit.framework.TestCase;

/**
 *         Date: 16/03/12
 *         Time: 10:12 AM
 */
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
    assertEquals("VcfBreakpoint: bar 100 foo 20 false false", b.toString());
    b = it.next();
    assertEquals("VcfBreakpoint: foo 20 bar 100 false false", b.toString());
    b = it.next();
    assertEquals("VcfBreakpoint: foo 22 bar 100 false false", b.toString());
    b = it.next();
    assertEquals("VcfBreakpoint: foo 22 foo 100 false false", b.toString());
    b = it.next();
    assertEquals("VcfBreakpoint: foo 22 foo 101 false false", b.toString());
    assertFalse(it.hasNext());
    b = bs.getMap().get("foo").get("bar").first();
    assertEquals("VcfBreakpoint: foo 20 bar 100 false false", b.toString());
    final List<String> names = bs.getChromosomes();
    assertEquals(2, names.size());
    assertEquals("bar", names.get(0));
    assertEquals("foo", names.get(1));
  }
}
