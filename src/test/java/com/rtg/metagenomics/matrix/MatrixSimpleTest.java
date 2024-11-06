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

import junit.framework.TestCase;

/**
 */
public class MatrixSimpleTest extends TestCase {

  public void test0() {
    final Matrix ma = new MatrixSimple(0);
    assertEquals("", ma.toString());
    final Matrix mb = new MatrixSimple(2);
    mb.set(0, 0, 0.5);
    mb.set(0, 1, 1.0);
    mb.set(1, 0, 0.2);
    mb.set(1, 1, 1.0);
    assertEquals(2, mb.rows());
    assertEquals(0.5, mb.get(0, 0));
    assertEquals(1.0, mb.get(0, 1));
    assertEquals(0.2, mb.get(1, 0));
    assertEquals(1.0, mb.get(1, 1));
  }

  public void test() {
    final Matrix ma = new MatrixSimple(3);
    assertTrue(ma.isSymmetric());
    final String exp = ""
      + "[0]  0.0000" + LS
      + "[1]  0.0000  0.0000" + LS
      + "[2]  0.0000  0.0000  0.0000" + LS
      ;
    assertEquals(exp, ma.toString());
    assertEquals(0.0, ma.get(0, 0));
    assertEquals(0.0, ma.get(0, 2));
    assertEquals(0.0, ma.get(2, 0));

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

    ma.set(2, 0, 2.5);
    final String exp2 = ""
      + "[0]  1.0000  0.0000  2.0000" + LS
      + "[1]  0.0000  0.0000  0.0000" + LS
      + "[2]  2.5000  0.0000  0.0000" + LS
      ;
    assertEquals(exp2, ma.toString());
    assertEquals(1.0, ma.get(0, 0));
    assertEquals(2.0, ma.get(0, 2));
    assertEquals(2.5, ma.get(2, 0));

    ma.incr(2, 0, 2.5);
    final String exp3 = ""
      + "[0]  1.0000  0.0000  2.0000" + LS
      + "[1]  0.0000  0.0000  0.0000" + LS
      + "[2]  5.0000  0.0000  0.0000" + LS
      ;
    assertEquals(exp3, ma.toString());
    assertEquals(1.0, ma.get(0, 0));
    assertEquals(2.0, ma.get(0, 2));
    assertEquals(5.0, ma.get(2, 0));
  }
}
