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
package com.rtg.mode;

import java.io.IOException;

import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class ScoringMatrixTest extends TestCase {

  private ScoringMatrix mMatrix;

  /**
   * Set up the matrix
   */
  @Override
  public void setUp() throws IOException, InvalidParamsException {
    mMatrix = new ProteinScoringMatrix("BLOSUM62");
  }

  /**
   * Tear down the matrix
   */
  @Override
  public void tearDown() {
    mMatrix = null;
  }


  public final void testScore() {
    check(4, Protein.A, Protein.A);
    check(-1, Protein.A, Protein.R);
    check(1, Protein.STOP, Protein.STOP);
    check(-1, Protein.K, Protein.D);
    check(2, Protein.M, Protein.L);
    check(0, Protein.X, Protein.A);
    check(-1, Protein.X, Protein.R);
    check(-2, Protein.X, Protein.C);
    check(-4, Protein.X, Protein.STOP);
  }

  private void check(final int sc, final Protein a, final Protein b) {
    assertEquals(sc, mMatrix.score(a.ordinal(), b.ordinal()));
    assertEquals(sc, mMatrix.score(b.ordinal(), a.ordinal()));
  }

  public final void testToStringStringBuilder() {
    final String str = mMatrix.toString();
    //System.err.println(str);
    assertTrue(str.startsWith("[-1, -4, 0, -1, -1, -1, -2,"));
    assertTrue(str.endsWith("-1, -2, -2, 0, -3, -1, 4]" + StringUtils.LS));
  }

}
