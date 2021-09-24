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
package com.rtg.metagenomics.metasnp;

import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class FrobeniusTest extends TestCase {

  public void testZero() {
    assertEquals(0.0, Frobenius.frobeniusNorm(new double[0][0]));
    assertEquals(0.0, Frobenius.frobeniusNorm(new double[1][1]));
    assertEquals(0.0, Frobenius.frobeniusNorm(new double[2][3]));
    assertEquals(0.0, Frobenius.frobeniusNorm(new double[3][2]));
  }

  public void testReal() {
    assertEquals(9.5393920141694564915262, Frobenius.frobeniusNorm(new double[][] {{1, 2, 3}, {4, 5, 6}}), 1e-10);
  }
  
  public void testDiff() {
    final double[][] a = {{1, 2, 3}, {4, 5, 6}};
    assertEquals(9.5393920141694564915262, Frobenius.frobeniusDistance(SimplePossibility.SINGLETON, new double[2][3], a), 1e-10);
    assertEquals(9.5393920141694564915262, Frobenius.frobeniusDistance(SimplePossibility.SINGLETON, a, new double[2][3]), 1e-10);
    assertEquals(0, Frobenius.frobeniusDistance(SimplePossibility.SINGLETON, a, a), 1e-10);
  }
}
