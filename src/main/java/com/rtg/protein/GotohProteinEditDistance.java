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
package com.rtg.protein;

import com.rtg.alignment.GotohEditDistance;
import com.rtg.mode.Protein;
import com.rtg.mode.ProteinScoringMatrix;

/**
 * This class subclasses the GotohEditDistance to work on protein.
 */
public class GotohProteinEditDistance extends GotohEditDistance {

  private final ProteinScoringMatrix mScoreMatrix;

  /**
   * Create a protein edit distance class with a particular scoring matrix.
   *
   * @param matrix one of the protein scoring matrices.
   */
  public GotohProteinEditDistance(final ProteinScoringMatrix matrix) {
    super((int) -(matrix.getGap() + matrix.getExpected()), (int) -matrix.getExtend(), 1, 0, false); //unknowns value should be ignored anyway (X in the matrix)
    mScoreMatrix = matrix;
  }

  @Override
  public boolean supportsEarlyTermination() {
    return false; //since protein can have negative penalties, early termination is problematic.
  }

  @Override
  protected int diagonalCost(final int templateResidue, final byte[] read, final int readPos) {
    final int readResidue = read[readPos - 1];
    // Off template matches are scored as X vs X according to the matrix
    return -mScoreMatrix.score(readResidue, templateResidue);
  }

  /**
   * Overridden to return amino acids.
   */
  @Override
  protected char residue(final byte[] a, final int pos) {
    return (char) Protein.pbase(a, pos);
  }

}
