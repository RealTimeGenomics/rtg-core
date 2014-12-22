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
public class AltFieldTest extends TestCase {

  public void testSomeMethod() {
    final String io = "##ALT=<ID=yo, Description=\"fun for the whole family\">";
    AltField f = new AltField(io);
    assertEquals(io.replaceAll(",\\s+", ","), f.toString());
    assertEquals("yo", f.getId());
    assertEquals("fun for the whole family", f.getDescription());
    assertEquals(f, f.superSet(f));
    final String io2 = "##ALT=<ID=%$&*,Description=\"oh rearry\">";
    final AltField f2 = new AltField(io2);
    assertEquals(io2, f2.toString());
    assertEquals("%$&*", f2.getId());
    assertEquals("oh rearry", f2.getDescription());
    assertEquals(f2, f2.superSet(f2));
    assertNull(f.superSet(f2));
    assertNull(f2.superSet(f));
  }
}
