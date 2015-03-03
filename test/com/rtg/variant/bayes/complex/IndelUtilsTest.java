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
package com.rtg.variant.bayes.complex;

import com.rtg.util.Utils;

import junit.framework.TestCase;

/**
 */
public class IndelUtilsTest extends TestCase {

  /**
   * Test method for {@link com.rtg.variant.bayes.complex.IndelUtils#flatten(double[], double)}.
   */
  public final void testFlatten1() {
    final double[] f = IndelUtils.flatten(new double[] {0.0, 0.5, 0.25, 0.25}, 0.3);
    assertEquals("[0.700, 0.150, 0.075, 0.075]", Utils.realFormat(f, 3));
  }

  /**
   * Test method for {@link com.rtg.variant.bayes.complex.IndelUtils#flatten(double[], double)}.
   */
  public final void testFlatten2() {
    final double[] f = IndelUtils.flatten(new double[] {0.0, 1.0}, 0.3);
    assertEquals("[0.700, 0.300]", Utils.realFormat(f, 3));
  }

  /**
   * Test method for {@link com.rtg.variant.bayes.complex.IndelUtils#unFlatten(double[])}.
   */
  public final void testUnFlatten1() {
    final double[] f = IndelUtils.unFlatten(new double[] {0.2, 0.2, 0.2, 0.2, 0.2});
    assertEquals("[0.000, 0.250, 0.250, 0.250, 0.250]", Utils.realFormat(f, 3));
  }

  /**
   * Test method for {@link com.rtg.variant.bayes.complex.IndelUtils#unFlatten(double[])}.
   */
  public final void testUnFlatten2() {
    final double[] f = IndelUtils.unFlatten(new double[] {0.5, 0.5});
    assertEquals("[0.000, 1.000]", Utils.realFormat(f, 3));
  }
}
