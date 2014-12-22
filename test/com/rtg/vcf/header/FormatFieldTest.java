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
package com.rtg.vcf.header;

import junit.framework.TestCase;

/**
 * Test class
 */
public class FormatFieldTest extends TestCase {

  public void testSomeMethod() {
    final String io = "##FORMAT=<ID=yo, Number=5, Type=Float, Description=\"fun for the whole family\">";
    final FormatField f = new FormatField(io);
    assertEquals(io.replaceAll(",\\s+", ","), f.toString());
    assertEquals("yo", f.getId());
    assertEquals(5, f.getNumber().getNumber());
    assertEquals(MetaType.FLOAT, f.getType());
    assertEquals("fun for the whole family", f.getDescription());
    assertEquals(f, f.superSet(f));
    final FormatField ioI = new FormatField(io.replaceAll("Float", "Integer"));
    assertEquals("yo", ioI.getId());
    assertEquals(5, ioI.getNumber().getNumber());
    assertEquals(MetaType.INTEGER, ioI.getType());
    assertEquals("fun for the whole family", ioI.getDescription());
    assertEquals(ioI, ioI.superSet(ioI));
    assertEquals(f, ioI.superSet(f));
    assertEquals(f, f.superSet(ioI));

    final String io2 = "##FORMAT=<ID=%$&*,Number=A,Type=Flag,Description=\"oh rearry\">";
    final FormatField f2 = new FormatField(io2);
    assertEquals(io2, f2.toString());
    assertEquals("%$&*", f2.getId());
    assertEquals(-1, f2.getNumber().getNumber());
    assertEquals(VcfNumberType.ALTS, f2.getNumber().getNumberType());
    assertEquals(MetaType.FLAG, f2.getType());
    assertEquals("oh rearry", f2.getDescription());
    assertEquals(f2, f2.superSet(f2));
    assertNull(f2.superSet(f));
    assertNull(f.superSet(f2));
  }
}
