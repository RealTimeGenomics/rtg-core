/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.cnv.preprocess;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class RegionColumnTest extends TestCase {

  public void test() {
    final RegionColumn col = new RegionColumn("col");
    assertEquals("col", col.getName());
    col.add("42:43+2");
    assertEquals(1, col.size());
    assertEquals("42:43-44", col.toString(0));
    col.add(0, "42:45+2");
    assertEquals("42:43-44", col.toString(1));
  }
}
