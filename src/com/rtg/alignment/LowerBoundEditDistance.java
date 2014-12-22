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
package com.rtg.alignment;

/**
 * Conservatively estimates the minimum alignment score and fails if it is more than a given maximum.
 *
 */
class LowerBoundEditDistance implements UnidirectionalEditDistance {

  private final LowerBoundEstimator mLBE;
  private final int[] mNoAlignmentPossible;

  /**
   * Default constructor
   * @param size the size of the <code>LowerBoundEstimator</code>
   * @param substitutionPenalty the penalty for a substitution
   * @param unknownsPenalty the penalty for an unknown nucleotide
   */
  protected LowerBoundEditDistance(int size, int substitutionPenalty, int unknownsPenalty) {
    mLBE = new LowerBoundEstimator(size, substitutionPenalty, unknownsPenalty);
    mNoAlignmentPossible = new int[ActionsHelper.ACTIONS_START_INDEX];
    mNoAlignmentPossible[ActionsHelper.ALIGNMENT_SCORE_INDEX] = Integer.MAX_VALUE;
  }

  /**
   * Lower bound estimator, using a 4-nt lock-step creeping count.
   * Returns null if the true alignment score might be less than or equal to <code>maxScore</code>.
   *
   * @param read the read
   * @param rlen read length
   * @param template the template
   * @param zeroBasedStart start position
   * @param maxScore maximum edit distance
   * @param maxShift the maximum allowed shift in start and end positions
   * @param cgLeft ignored
   * @return null, or an alignment score of <code>Integer.MAX_VALUE</code>.
   */
  @Override
  public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {

    final int lb = mLBE.calcLB(read, rlen, template, zeroBasedStart, maxScore, maxShift);
    if (lb > maxScore) {
      mNoAlignmentPossible[ActionsHelper.TEMPLATE_START_INDEX] = zeroBasedStart;
      return mNoAlignmentPossible;
    }

    return null;
  }

  @Override
  public void logStats() {
  }

  @Override
  public int[] calculateEditDistanceFixedBoth(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos,
      int templateEndPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int[] calculateEditDistanceFixedEnd(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateExpectedStartPos,
      int templateEndPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int[] calculateEditDistanceFixedStart(byte[] read, int readStartPos, int readEndPos, byte[] template, int templateStartPos, int maxScore, int maxShift) {
    throw new UnsupportedOperationException();
  }
}
