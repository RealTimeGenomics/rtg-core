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
package com.rtg.variant.cnv.region;

import junit.framework.TestCase;

/**
 */
public class EmptyRegionTest extends TestCase {

  public void test() {
    final Region r = EmptyRegion.EMPTY_REGION;
    assertFalse(r.contains(0));
    assertFalse(r.contains(Integer.MAX_VALUE));
    assertFalse(r.contains(Integer.MIN_VALUE));
    assertEquals("EmptyRegion", r.toString());
  }
}
