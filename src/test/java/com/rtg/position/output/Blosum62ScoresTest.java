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
import junit.framework.TestCase;


/**
 */
public class Blosum62ScoresTest extends TestCase {

  private static final double NEG = Double.NEGATIVE_INFINITY;

  public void test() {
    final Blosum62Scores gp = new Blosum62Scores(25);
    gp.globalIntegrity();
    //System.err.println(gp.toString());
    assertEquals(-15, gp.minDelta());
    assertEquals(+15, gp.maxDelta());

    assertEquals(-0.0, gp.score(0, 0, 0, 0));

    assertEquals(-11.0, gp.score(0, 0, 0, 1)); //0 - -1 -1 = 0, 1 --1 -1 = 1
    assertEquals(-11.0, gp.score(0, 1, 0, 0));

    assertEquals(-7.8135, gp.score(0, 15, 0, 15));
    assertEquals(-18.2926, gp.score(0, 14, 0, 15));

    assertEquals(-8.3344, gp.score(0, 16, 0, 16));
    assertEquals(NEG, gp.score(0, 16, 0, 15));
    assertEquals(NEG, gp.score(0, 15, 0, 16));

    assertEquals(-13.0225, gp.score(0, 25, 0, 25));
    assertEquals(NEG, gp.score(0, 24, 0, 25));
    assertEquals(NEG, gp.score(0, 25, 0, 24));

    assertEquals(NEG, gp.score(0, 26, 0, 26));

    assertEquals(-0.0, gp.scoreMax(0, 0, -1, -1));
    assertEquals(-7.8135, gp.scoreMax(0, 15, -1, -1));
    assertEquals(-13.0225, gp.scoreMax(0, 25, -1, -1));
    assertEquals(NEG, gp.scoreMax(0, 26, -1, -1));

    assertEquals("Blosum62Probabilities", gp.toString());
  }
}
