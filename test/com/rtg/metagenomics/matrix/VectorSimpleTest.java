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


import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 */
public class VectorSimpleTest extends TestCase {

  public void test0() {
    final VectorSimple ve = new VectorSimple(0);
    ve.globalIntegrity();
    Exam.assertEquals("", ve.toString());
    assertEquals(0, ve.dimension());
  }

  public void test() {
    final VectorSimple ve = new VectorSimple(3);
    ve.globalIntegrity();
    Exam.assertEquals("  0.0000  0.0000  0.0000", ve.toString());
    Exam.assertEquals(0.0, ve.get(0));
    Exam.assertEquals(0.0, ve.get(1));
    Exam.assertEquals(0.0, ve.get(2));

    ve.set(0, 0.5);
    ve.set(1, 1.0);
    ve.set(2, 2.0);
    Exam.assertEquals("  0.5000  1.0000  2.0000", ve.toString());
    Exam.assertEquals(0.5, ve.get(0));
    Exam.assertEquals(1.0, ve.get(1));
    Exam.assertEquals(2.0, ve.get(2));

    ve.set(2, 2.5);
    Exam.assertEquals("  0.5000  1.0000  2.5000", ve.toString());
    Exam.assertEquals(0.5, ve.get(0));
    Exam.assertEquals(1.0, ve.get(1));
    Exam.assertEquals(2.5, ve.get(2));

    ve.incr(2, 2.5);
    Exam.assertEquals("  0.5000  1.0000  5.0000", ve.toString());
    Exam.assertEquals(0.5, ve.get(0));
    Exam.assertEquals(1.0, ve.get(1));
    Exam.assertEquals(5.0, ve.get(2));
  }
}
