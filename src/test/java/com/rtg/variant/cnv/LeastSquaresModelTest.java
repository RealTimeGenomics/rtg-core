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
package com.rtg.variant.cnv;

import static com.rtg.util.StringUtils.LS;

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import junit.framework.TestCase;

/**
 */
public class LeastSquaresModelTest extends TestCase {

  public void test() throws IOException {
    final double[] ratios = {0.5, 0.75, 0.25, 0.4, 0.6, 0.3, 0.7, 1.0, 1.2, 0.8, 0.7, 1.3, 0.6, 1.4};
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final OutputStream ps = new BufferedOutputStream(bos);
    final GnuOutput gnu = new GnuOutput(ps);
    final LeastSquaresModel lsm = new LeastSquaresModel("test", 2, 1, 0.5, false, gnu);
    for (final double ratio : ratios) {
      lsm.add(ratio);
    }
    lsm.globalIntegrity();
    lsm.fit();
    ps.close();
    final String exp = ""
      + "-0.25 0.750" + LS
      + "13.25 0.750" + LS
      + "" + LS
      + "-0.25 0.391" + LS
      + "-0.25 1.109" + LS
      + "" + LS
      + "13.25 0.391" + LS
      + "13.25 1.109" + LS
      + "" + LS
      ;
    assertEquals(exp, bos.toString());
  }

  public void testZero() throws IOException {
    final double[] ratios = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final OutputStream ps = new BufferedOutputStream(bos);
    final GnuOutput gnu = new GnuOutput(ps);
    final LeastSquaresModel lsm = new LeastSquaresModel("test", 2, 1, 0.0, false, gnu);
    for (final double ratio : ratios) {
      lsm.add(ratio);
    }
    lsm.globalIntegrity();
    lsm.fit();
    ps.close();
    final String exp = ""
      + "-0.25 0.000" + LS
      + "13.25 0.000" + LS
      + "" + LS
      + "-0.25 0.000" + LS
      + "-0.25 0.000" + LS
      + "" + LS
      + "13.25 0.000" + LS
      + "13.25 0.000" + LS
      + "" + LS
      ;
    assertEquals(exp, bos.toString());
  }

  public void testZeroThousand() throws IOException {
    final double[] ratios = new double[1000];
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final OutputStream ps = new BufferedOutputStream(bos);
    final GnuOutput gnu = new GnuOutput(ps);
    final LeastSquaresModel lsm = new LeastSquaresModel("test", 2, 1, 0.0, false, gnu);
    for (final double ratio : ratios) {
      lsm.add(ratio);
    }
    lsm.globalIntegrity();
    lsm.fit();
    ps.close();
    final String exp = ""
      + "-0.25 0.000" + LS
      + "999.25 0.000" + LS
      + "" + LS
      + "-0.25 0.000" + LS
      + "-0.25 0.000" + LS
      + "" + LS
      + "999.25 0.000" + LS
      + "999.25 0.000" + LS
      + "" + LS
      ;
    assertEquals(exp, bos.toString());
  }

  public void testZeroInfinity() throws IOException {
    final double[] ratios = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final OutputStream ps = new BufferedOutputStream(bos);
    final GnuOutput gnu = new GnuOutput(ps);
    final LeastSquaresModel lsm = new LeastSquaresModel("test", 2, 1, 0.0, false, gnu);
    for (final double ratio : ratios) {
      lsm.add(ratio);
    }
    lsm.globalIntegrity();
    lsm.fit();
    ps.close();
    final String exp = ""
      + "-0.25 1.000" + LS
      + "5.25 1.000" + LS
      + "" + LS
      + "-0.25 1.000" + LS
      + "-0.25 1.000" + LS
      + "" + LS
      + "5.25 1.000" + LS
      + "5.25 1.000" + LS
      + "" + LS
      + "5.75 0.000" + LS
      + "13.25 0.000" + LS
      + "" + LS
      + "5.75 0.000" + LS
      + "5.75 0.000" + LS
      + "" + LS
      + "13.25 0.000" + LS
      + "13.25 0.000" + LS
      + "" + LS
      ;
    assertEquals(exp, bos.toString());
  }
}
