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
    mNoAlignmentPossible[ActionsHelper.TEMPLATE_START_INDEX] = Integer.MAX_VALUE;
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
