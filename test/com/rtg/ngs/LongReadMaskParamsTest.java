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

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class LongReadMaskParamsTest extends TestCase {

  public void testGetters() {
    final int readLength = 35;
    final LongReadMaskParams lrmp = new LongReadMaskParams(25, 1, 2, 3);
    assertEquals(25, lrmp.getWordSize());
    assertEquals(1, lrmp.getSubstitutions());
    assertEquals(2, lrmp.getIndels());
    assertEquals(3, lrmp.getIndelLength());

    assertTrue(lrmp.closed());
    assertTrue(lrmp.isValid(readLength));
    assertEquals("Long read Mask: w=25 s=1 i=2 l=3", lrmp.toString());

    try {
      lrmp.maskFactory(readLength);
      fail();
    } catch (UnsupportedOperationException ex) {
      assertEquals("Not supported, ever", ex.getMessage());
    }
  }


}
