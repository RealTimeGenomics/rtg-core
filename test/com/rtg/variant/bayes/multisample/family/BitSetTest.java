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

package com.rtg.variant.bayes.multisample.family;

import junit.framework.TestCase;

/**
 */
public class BitSetTest extends TestCase {

  public void test() {
    final BitSet bs = BitSet.DNA_SET;
    bs.integrity();
    assertEquals(4, bs.length());
    assertEquals("A:C:G:T", bs.toString());
    assertEquals("empty", bs.toString(0));
    assertEquals("A", bs.toString(1 << 0));
    assertEquals("C", bs.toString(1 << 1));
    assertEquals("G", bs.toString(1 << 2));
    assertEquals("T", bs.toString(1 << 3));

    assertEquals(0, bs.toSet());
    assertEquals("A", bs.toString(bs.toSet(0)));
    assertEquals("C", bs.toString(bs.toSet(1)));
    assertEquals("G", bs.toString(bs.toSet(2)));
    assertEquals("T", bs.toString(bs.toSet(3)));

    assertTrue(bs.contains(1, 0));
    assertFalse(bs.contains(1, 1));
    assertFalse(bs.contains(1, 2));
    assertFalse(bs.contains(1, 3));

    assertTrue(bs.contains(15, 0));
    assertTrue(bs.contains(15, 1));
    assertTrue(bs.contains(15, 2));
    assertTrue(bs.contains(15, 3));

    assertEquals(12, bs.complement(3));
    assertEquals(15, bs.complement(0));
  }

  public void test1() {
    final BitSet bs = new BitSet("A");
    bs.integrity();
    assertEquals(1, bs.length());
    assertEquals("A", bs.toString());
    assertEquals("empty", bs.toString(0));
    assertEquals("A", bs.toString(1 << 0));

    assertEquals(0, bs.toSet());
    assertEquals("A", bs.toString(bs.toSet(0)));

    assertEquals(0, bs.complement(1));
    assertEquals(1, bs.complement(0));
  }

  //32 entries - special case for bit handling
  public void test32() {
    final String[] s = new String[32];
    for (int i = 0; i < 32; i++) {
      s[i] = "" + i;
    }
    final BitSet bs = new BitSet(s);
    bs.integrity();
    assertEquals(32, bs.length());
    assertEquals("empty", bs.toString(0));
    for (int i = 0; i < 32; i++) {
      assertEquals(s[i], bs.toString(1 << i));
    }

    assertEquals(s[0] + ":" + s[31], bs.toString(bs.toSet(0, 31)));
    assertEquals(0, bs.toSet());

    assertEquals(Integer.MAX_VALUE, bs.complement(Integer.MIN_VALUE));
    assertEquals(-1, bs.complement(0));

  }
}
