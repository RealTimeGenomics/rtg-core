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
