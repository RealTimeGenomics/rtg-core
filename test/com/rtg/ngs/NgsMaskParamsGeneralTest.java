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
package com.rtg.ngs;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class NgsMaskParamsGeneralTest extends TestCase {

  NgsMaskParamsGeneral getParams(final int wordsize, final int substitutions, final int indels, final int indelLength) {
    return new NgsMaskParamsGeneral(wordsize, substitutions, indels, indelLength, false);
  }

  public void testEquals() {
    final int readLength = 36;
    final NgsMaskParams a1 = getParams(12, 2, 1, 1);
    assertTrue(a1.isValid(readLength));
    final NgsMaskParams a2 = getParams(12, 2, 1, 1);
    assertTrue(a2.isValid(readLength));
    final NgsMaskParams z = getParams(12, 2, 2, 2);
    assertTrue(z.isValid(readLength));
    final NgsMaskParams b = getParams(12, 2, 2, 1);
    assertTrue(b.isValid(readLength));
    final NgsMaskParams c = getParams(12, 3, 1, 1);
    assertTrue(c.isValid(readLength));
    final NgsMaskParams d = getParams(13, 2, 1, 1);
    assertTrue(d.isValid(readLength));
    TestUtils.equalsHashTest(new NgsMaskParams[][] {{a1, a2}, {z}, {b}, {c}, {d}});
  }


  public void test() throws Exception {
    final NgsMaskParamsGeneral a = getParams(12, 3, 2, 1);
    a.integrity();
    assertEquals("General Mask: w=12 s=3 i=2 l=1", a.toString());
    assertNotNull(a.maskFactory(36));
    assertTrue(a.closed());
    a.close();
    assertTrue(a.closed());
  }

  public void testInvalid() {
    final NgsMaskParams a = getParams(32, 1, 1, 1);
    assertFalse(a.isValid(33));
  }
}
