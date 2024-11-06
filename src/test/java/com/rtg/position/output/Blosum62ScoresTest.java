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
