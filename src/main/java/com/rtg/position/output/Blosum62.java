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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import com.rtg.mode.Protein;
import com.rtg.mode.ScoringMatrix;
import com.rtg.util.Resources;
import com.rtg.util.diagnostic.SlimException;

/**
 * Holds constants from the Blosum62 matrix using gapped scoring.
 * These are collected together here for convenience and to ensure consistent values are being used.
 * <pre>
 * E = K * m * n * exp (- LAMBDA * raw score)
 * </pre>
 */
public final class Blosum62 extends ScoringMatrix {

  /** K in expect equation. */
  public static final double K = 0.0410;

  /** LAMBDA in expect equation. */
  public static final double LAMBDA = 0.267;

  /** Average score for two identical amino acids. */
  public static final double HIT = 4.5;

  /** Average score for two different amino acids. */
  public static final double MISS = -1.0;

  /** Score on opening a gap. */
  public static final double GAP = -10.0;

  /** Score on extending a gap. */
  public static final double EXTEND = -1.0;

  /** Expected value for score. */
  public static final double EXPECTED = -0.5209;

  /** Single instance. */
  public static final Blosum62 SINGLETON = new Blosum62();

  /**
   * Blosum62 matrix for protein distance.  Use the SINGLETON to get an instance.
   */
  public Blosum62() {
    final int len = Protein.values().length;
    mScores = new int[len][len];
    try {
      try (BufferedReader re = new BufferedReader(new InputStreamReader(Resources.getResourceAsStream("com/rtg/mode/BLOSUM62")))) {
        parse(re);
      }
    } catch (final IOException e) {
      throw new SlimException(e);
    }
    assert integrity();
  }

}
