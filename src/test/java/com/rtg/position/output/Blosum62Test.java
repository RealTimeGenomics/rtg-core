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
package com.rtg.position.output;

import static com.rtg.mode.Protein.A;
import static com.rtg.mode.Protein.C;
import static com.rtg.mode.Protein.D;
import static com.rtg.mode.Protein.K;
import static com.rtg.mode.Protein.L;
import static com.rtg.mode.Protein.M;
import static com.rtg.mode.Protein.R;
import static com.rtg.mode.Protein.STOP;
import static com.rtg.mode.Protein.X;

import com.rtg.mode.Protein;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class Blosum62Test extends TestCase {

  /**
   * Test method for {@link Blosum62#score(int, int)}.
   */
  public final void testScore() {
    check(4, A, A);
    check(-1, A, R);
    check(1, STOP, STOP);
    check(-1, K, D);
    check(2, M, L);
    check(0, X, A);
    check(-1, X, R);
    check(-2, X, C);
    check(-4, X, STOP);
  }

  private void check(final int sc, final Protein a, final Protein b) {
    assertEquals(sc, Blosum62.SINGLETON.score(a.ordinal(), b.ordinal()));
    assertEquals(sc, Blosum62.SINGLETON.score(b.ordinal(), a.ordinal()));
  }

  public final void testToStringStringBuilder() {
    final String str = Blosum62.SINGLETON.toString();
    //System.err.println(str);
    assertTrue(str.startsWith("[-1, -4, 0, -1, -1, -1, -2,"));
    assertTrue(str.endsWith("-1, -2, -2, 0, -3, -1, 4]" + StringUtils.LS));
  }

}
