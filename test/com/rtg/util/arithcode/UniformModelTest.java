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

package com.rtg.util.arithcode;

import junit.framework.TestCase;

/**
 */
public class UniformModelTest extends TestCase {

  public void test() {
    final UniformModel mo = new UniformModel(3);
    assertEquals(3, mo.totalCount());

    assertEquals(0, mo.pointToSymbol(0));
    assertEquals(1, mo.pointToSymbol(1));
    assertEquals(2, mo.pointToSymbol(2));

    //assertEquals(3, mo.pointToSymbol(3));
    final int[] result = new int[3];
    mo.interval(0, result);
    assertEquals(0, result[0]);
    assertEquals(1, result[1]);
    assertEquals(3, result[2]);

    mo.interval(1, result);
    assertEquals(1, result[0]);
    assertEquals(2, result[1]);
    assertEquals(3, result[2]);

    mo.interval(2, result);
    assertEquals(2, result[0]);
    assertEquals(3, result[1]);
    assertEquals(3, result[2]);

    assertEquals(256, UniformModel.MODEL.totalCount());
  }
}
