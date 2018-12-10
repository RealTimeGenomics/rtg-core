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
public class IntColumnTest extends TestCase {

  public void test() {
    final IntColumn col = new IntColumn("col");
    assertEquals("col", col.getName());
    col.add(42);
    assertEquals(1, col.size());
    assertEquals(42.0, col.get(0));
    assertEquals("42", col.toString(0));
    assertEquals(42.0, col.mean());
    assertEquals(42.0, col.median());
    col.remove(0);
    assertEquals(0, col.size());
    col.add(1);
    col.add(2);
    col.add(3);
    col.add(4);
    assertEquals(4, col.size());
    col.remove(1);  // 1,3,4
    assertEquals(3, col.size());
    assertEquals(1, (int) col.get(0));
    assertEquals(3, (int) col.get(1));
    col.add(0, 3);  // 3,1,3,4
    assertEquals(3, (int) col.get(0));
    assertEquals(1, (int) col.get(1));
    assertEquals(3, (int) col.get(2));
    assertEquals(4, (int) col.get(3));
    col.add(3, 5);  // 3,1,3,5,4
    assertEquals(5, col.size());
    assertEquals(4, (int) col.get(4));
  }
}
