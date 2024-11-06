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
package com.rtg.alignment;

import com.rtg.mode.Protein;
import com.rtg.mode.ProteinScoringMatrix;

/**
 * Provides validation methods for the packed actions arrays returned by edit distance classes.
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
