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
public class NgsMaskParamsExplicitTest extends TestCase {

  NgsMaskParamsExplicit getParams(final String mask) {
    return new NgsMaskParamsExplicit(mask);
  }

  public void testEquals() throws Exception {

    final NgsMaskParams a1 = getParams("SplitL4w4s0e0");
    final NgsMaskParams a2 = getParams("SplitL4w4s0e0");
    final NgsMaskParams b = getParams("SplitL4w2s1e1b");
    TestUtils.equalsHashTest(new NgsMaskParams[][] {{a1, a2}, {b}});
    a1.close();
    a2.close();
    b.close();
  }


  public void test() throws Exception {
    final NgsMaskParamsExplicit a = getParams("SplitL4w4s0e0");
    a.integrity();
    assertEquals("Mask:SplitL4w4s0e0", a.toString());
    assertNotNull(a.maskFactory(33));
    assertTrue(a.closed());
    a.close();
    assertTrue(a.closed());
  }
}
