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

import junit.framework.TestCase;

/**
 */
public class ProbAlphaSimpleBetaTest extends TestCase {
  public void test() {
    ProbAlphaSimpleBeta beta = new ProbAlphaSimpleBeta(new double[] {0.01, 0.0002, 0.0003});
    assertEquals(0.01 * 0.0002 * 0.0003, beta.pAlpha(0, new int[] {1, 1, 1}), 1e-8);
    assertEquals(0.99 * 0.0002 * 0.0003, beta.pAlpha(0, new int[] {0, 1, 1}), 1e-8);
    assertEquals(0.01 * 0.0002 * 0.9997, beta.pAlpha(0, new int[] {1, 1, 0}), 1e-8);
  }
}
