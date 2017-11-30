/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.alignment;

import com.rtg.mode.Protein;
import com.rtg.mode.ProteinScoringMatrix;

/**
 * Provides validation methods for the packed actions arrays returned by EditDistance classes.
 *
 */
public class ProteinActionsValidator extends ActionsValidator {

  private final ProteinScoringMatrix mScoreMatrix;

  /**
   * Construct a actions validator for Protein edit distance classes.
   *
   * @param matrix the protein scoring matrix.
   */
  ProteinActionsValidator(final ProteinScoringMatrix matrix) {
    super((int) -matrix.getGap(), (int) -matrix.getExtend(), null, null, (byte) Protein.X.ordinal());
    mScoreMatrix = matrix;
  }

  @Override
  public boolean isValid(int[] actions, byte[] read, int rlen, byte[] template, boolean rc, int maxScore) {
    assert !rc;
    return isValid(actions, read, rlen, template, maxScore);
  }

  @Override
  protected int scoreSame(final byte base, final byte tbase) {
    return -mScoreMatrix.score(base, tbase);
  }

  @Override
  protected int scoreSubs(final byte base, final byte tbase) {
    return -mScoreMatrix.score(base, tbase);
  }
}
