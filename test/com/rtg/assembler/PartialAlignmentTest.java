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

package com.rtg.assembler;

import junit.framework.TestCase;

/**
 */
public class PartialAlignmentTest extends TestCase {
  public void test() {
    PartialAlignment pa = new PartialAlignment(2, 10, 20, 30, 40, 50);
    assertFalse(pa.equals(null));
    assertFalse(pa.equals("ASDF"));
    assertTrue(pa.equals(pa));
    PartialAlignment pa2 = new PartialAlignment(3, 10, 20, 30, 40, 50);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(3, 11, 20, 30, 40, 50);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(2, 10, 22, 30, 40, 50);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(2, 10, 20, 32, 40, 50);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(2, 10, 20, 30, 42, 50);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(2, 10, 20, 30, 40, 52);
    assertFalse(pa.equals(pa2));
    pa2 = new PartialAlignment(2, 10, 20, 30, 40, 50);
    assertTrue(pa.equals(pa2));
    assertEquals(pa.hashCode(), pa2.hashCode());
    assertEquals(-916452864, pa.hashCode());
  }
  public void testGetters() {
    PartialAlignment pa = new PartialAlignment(2, 10, 20, 30, 40, 50);
    assertEquals(2, pa.getAlignmentScore());
    assertEquals(10, pa.getReadStart());
    assertEquals(20, pa.getReadEnd());
    assertEquals(30, pa.getContig());
    assertEquals(40, pa.getContigStart());
    assertEquals(50, pa.getContigEnd());
    assertEquals("contig=30 [40,50] read [10,20] score=2", pa.toString());
  }
}
