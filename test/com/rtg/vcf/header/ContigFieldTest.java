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
public class ContigFieldTest extends TestCase {

  public void testSomeMethod() {
    final String io = "##contig=<ID=yo,length=20>";
    final ContigField f = new ContigField(io);
    assertEquals(io.replaceAll(",\\s+", ","), f.toString());
    assertEquals("yo", f.getId());
    assertEquals(Integer.valueOf(20), f.getLength());
    assertEquals(f, f.superSet(f));
    final String io2 = "##contig=<ID=%$&*,length=17>";
    final ContigField f2 = new ContigField(io2);
    assertEquals(io2, f2.toString());
    assertEquals("%$&*", f2.getId());
    assertEquals(Integer.valueOf(17), f2.getLength());
    assertEquals(f2, f2.superSet(f2));
    assertNull(f.superSet(f2));
    assertNull(f2.superSet(f));

    final String io3 = "##contig=<ID=yo,assembly=b37>";
    final ContigField f3 = new ContigField(io3);
    assertEquals(io3, f3.toString());
    assertEquals("yo", f3.getId());
    assertEquals(null, f3.getLength());
    assertEquals(f3, f3.superSet(f3));
    assertNull(f2.superSet(f3));
    final ContigField f4 = f.superSet(f3);
    assertEquals("yo", f4.getId());
    assertEquals(Integer.valueOf(20), f4.getLength());
    assertEquals("##contig=<ID=yo,length=20,assembly=b37>", f4.toString());
    final ContigField f5 = f3.superSet(f);
    assertEquals(f4, f5);
  }
}
