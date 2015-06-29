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
package com.rtg.util;

import junit.framework.TestCase;

/**
 */
public class BoundedDoubleTest extends TestCase {

  public void testConstruction() {
    final BoundedDouble bd = new BoundedDouble(0.45, 0.2, 0.63);
    assertEquals(0.2, bd.getLow());
    assertEquals(0.45, bd.getValue());
    assertEquals(0.63, bd.getHigh());
    assertEquals("0.45", bd.toString());
  }
}
