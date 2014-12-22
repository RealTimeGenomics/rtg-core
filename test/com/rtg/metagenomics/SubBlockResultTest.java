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
package com.rtg.metagenomics;

import com.rtg.metagenomics.matrix.VectorSimple;

import junit.framework.TestCase;

/**
 */
public class SubBlockResultTest extends TestCase {


  public void test() {
    final VectorSimple x = new VectorSimple(3);
    final VectorSimple varianceLog = new VectorSimple(3);
    final VectorSimple likelihood = new VectorSimple(3);
    final SubBlockResult res = new SubBlockResult(x, varianceLog,  likelihood, 1.5);
    assertEquals(3, res.getR().dimension());
    assertEquals(3, res.getVarianceLog().dimension());
    assertEquals(3, res.getLikelihoods().dimension());
    assertEquals(1.5, res.getL());
  }
}
