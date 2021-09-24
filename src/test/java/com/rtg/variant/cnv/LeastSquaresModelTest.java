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
