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
package com.rtg.protein;

import com.rtg.alignment.ActionsHelper;
import com.rtg.alignment.EditDistance;
import com.rtg.alignment.GotohEditDistance;
import com.rtg.mode.ProteinScoringMatrix;

/**
 * This class subclasses the GotohEditDistance to work on protein.
 */
public class GotohProteinEditDistance extends GotohEditDistance implements EditDistance {

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
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, boolean rc, int maxScore, int maxShift, boolean cgLeft) {
    if (rc) {
      throw new RuntimeException("Gotoh does not support reverse complement (ProteinEditDistance)");
    }
    return calculateEditDistance(read, rlen, template, zeroBasedStart, maxScore, maxShift, cgLeft);
  }

//  @Override
//  public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, boolean rc,
//      int maxScore, int maxShift) {
//    if (rc) {
//      throw new RuntimeException("Gotoh does not support reverse complement (ProteinEditDistance)");
//    }
//    return calculateEditDistanceFixedStart(read, readStartPos, readEndPos, template, templateStartPos, maxScore, maxShift);
//  }

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
    return (char) ActionsHelper.pbase(a, pos);
  }

}
