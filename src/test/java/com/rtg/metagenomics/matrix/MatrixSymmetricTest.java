/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
