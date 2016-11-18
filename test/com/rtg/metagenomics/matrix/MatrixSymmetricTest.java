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

import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 */
public class MatrixSymmetricTest extends TestCase {

  public void test0() {
    final Matrix ma = new MatrixSymmetric(0);
    Exam.assertEquals("", ma.toString());
  }

  public void test() {
    final Matrix ma = new MatrixSymmetric(3);
    Exam.assertTrue(ma.isSymmetric());
    final String exp = ""
      + "[0]  0.0000" + LS
      + "[1]  0.0000  0.0000" + LS
      + "[2]  0.0000  0.0000  0.0000" + LS
      ;
    Exam.assertEquals(exp, ma.toString());
    Exam.assertEquals(0.0, ma.get(0, 0));
    Exam.assertEquals(0.0, ma.get(0, 2));
    Exam.assertEquals(0.0, ma.get(2, 0));

    ma.set(0, 0, 1.0);
    ma.set(0, 2, 2.0);
    final String exp1 = ""
      + "[0]  1.0000" + LS
      + "[1]  0.0000  0.0000" + LS
      + "[2]  2.0000  0.0000  0.0000" + LS
      ;
    Exam.assertEquals(exp1, ma.toString());
    Exam.assertEquals(1.0, ma.get(0, 0));
    Exam.assertEquals(2.0, ma.get(0, 2));
    Exam.assertEquals(2.0, ma.get(2, 0));

    ma.set(2, 0, 2.5);
    final String exp2 = ""
      + "[0]  1.0000" + LS
      + "[1]  0.0000  0.0000" + LS
      + "[2]  2.5000  0.0000  0.0000" + LS
      ;
    Exam.assertEquals(exp2, ma.toString());
    Exam.assertEquals(1.0, ma.get(0, 0));
    Exam.assertEquals(2.5, ma.get(0, 2));
    Exam.assertEquals(2.5, ma.get(2, 0));

    ma.incr(2, 0, 2.5);
    final String exp3 = ""
      + "[0]  1.0000" + LS
      + "[1]  0.0000  0.0000" + LS
      + "[2]  5.0000  0.0000  0.0000" + LS
      ;
    Exam.assertEquals(exp3, ma.toString());
    Exam.assertEquals(1.0, ma.get(0, 0));
    Exam.assertEquals(5.0, ma.get(0, 2));
    Exam.assertEquals(5.0, ma.get(2, 0));

    ma.incr(0, 2, 2.5);
    final String exp4 = ""
      + "[0]  1.0000" + LS
      + "[1]  0.0000  0.0000" + LS
      + "[2]  7.5000  0.0000  0.0000" + LS
      ;
    Exam.assertEquals(exp4, ma.toString());
    Exam.assertEquals(1.0, ma.get(0, 0));
    Exam.assertEquals(7.5, ma.get(0, 2));
    Exam.assertEquals(7.5, ma.get(2, 0));
}
}
