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
package com.rtg.metagenomics.matrix;

import static com.rtg.util.StringUtils.LS;

import junit.framework.TestCase;

/**
 */
public class MatrixTransposeTest extends TestCase {

  public void test0() {
    final Matrix ma = new MatrixTranspose(new MatrixSimple(0));
    assertEquals("", ma.toString());
  }

  public void test() {
    final Matrix ma = new MatrixSimple(3);
    assertTrue(ma.isSymmetric());
    assertEquals(3, ma.size());

    ma.set(0, 0, 1.0);
    ma.set(0, 2, 2.0);
    final String exp1 = ""
      + "[0]  1.0000  0.0000  2.0000" + LS
      + "[1]  0.0000  0.0000  0.0000" + LS
      + "[2]  0.0000  0.0000  0.0000" + LS
      ;
    assertEquals(exp1, ma.toString());
    assertEquals(1.0, ma.get(0, 0));
    assertEquals(2.0, ma.get(0, 2));
    assertEquals(0.0, ma.get(2, 0));

    final Matrix mt = new MatrixTranspose(ma);
    assertFalse(mt.isSymmetric());
    assertEquals(3, mt.size());

    final String exp2 = ""
      + "[0]  1.0000  0.0000  0.0000" + LS
      + "[1]  0.0000  0.0000  0.0000" + LS
      + "[2]  2.0000  0.0000  0.0000" + LS
      ;
    assertEquals(exp2, mt.toString());
    assertEquals(1.0, mt.get(0, 0));
    assertEquals(0.0, mt.get(0, 2));
    assertEquals(2.0, mt.get(2, 0));

    mt.incr(0, 2, 2.5);
    mt.set(1, 0, 4.2);
    final String exp3 = ""
      + "[0]  1.0000  0.0000  2.5000" + LS
      + "[1]  4.2000  0.0000  0.0000" + LS
      + "[2]  2.0000  0.0000  0.0000" + LS
      ;
    assertEquals(exp3, mt.toString());
    assertEquals(1.0, mt.get(0, 0));
    assertEquals(4.2, mt.get(1, 0));
    assertEquals(2.0, mt.get(2, 0));
    assertEquals(2.5, mt.get(0, 2));
  }
}
