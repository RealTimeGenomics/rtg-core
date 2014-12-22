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

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 * Test class
 */
public class VcfNumberTest extends TestCase {

  public void testSomeMethod() {
    final VcfNumber a = new VcfNumber("A");
    assertEquals("A", a.toString());
    assertEquals(VcfNumberType.ALTS, a.getNumberType());
    assertEquals(-1, a.getNumber());
    final VcfNumber a2 = new VcfNumber("G");
    assertEquals("G", a2.toString());
    assertEquals(VcfNumberType.GENOTYPES, a2.getNumberType());
    assertEquals(-1, a2.getNumber());
    final VcfNumber a3 = new VcfNumber(".");
    assertEquals(".", a3.toString());
    assertEquals(VcfNumberType.UNKNOWN, a3.getNumberType());
    assertEquals(-1, a3.getNumber());
    final VcfNumber a4 = new VcfNumber("13");
    assertEquals("13", a4.toString());
    assertEquals(VcfNumberType.INTEGER, a4.getNumberType());
    assertEquals(13, a4.getNumber());
    final VcfNumber a5 = new VcfNumber("0");
    assertEquals("0", a5.toString());
    assertEquals(VcfNumberType.INTEGER, a5.getNumberType());
    assertEquals(0, a5.getNumber());
    final VcfNumber a6 = new VcfNumber("R");
    assertEquals("R", a6.toString());
    assertEquals(VcfNumberType.REF_ALTS, a6.getNumberType());
    assertEquals(-1, a6.getNumber());

    assertTrue(a.equals(a));
    assertTrue(a2.equals(a2));
    assertTrue(a3.equals(a3));
    assertTrue(a4.equals(a4));
    assertTrue(a5.equals(a5));
    assertTrue(a6.equals(a6));
    assertFalse(a4.equals(a5));
    assertFalse(a.equals(a5));
    assertFalse(a.equals(new VcfNumber("-1")));
    assertTrue(a4.equals(new VcfNumber("13")));
    assertTrue(a5.equals(new VcfNumber("0")));
  }

  public void testHashCode() {
    TestUtils.equalsHashTest(new VcfNumber[][] {new VcfNumber[] {new VcfNumber("A"), new VcfNumber("A")},
                                                new VcfNumber[] {new VcfNumber("R"), new VcfNumber("R")},
                                                new VcfNumber[] {new VcfNumber("G"), new VcfNumber("G")},
                                                new VcfNumber[] {new VcfNumber("."), new VcfNumber(".")},
                                                new VcfNumber[] {new VcfNumber("13"), new VcfNumber("13")},
                                                new VcfNumber[] {new VcfNumber("0"), new VcfNumber("0")},
                                                new VcfNumber[] {new VcfNumber("-1"), new VcfNumber("-1")},
    });
  }
}
