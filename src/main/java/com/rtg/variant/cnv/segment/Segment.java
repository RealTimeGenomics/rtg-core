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
package com.rtg.variant.cnv.segment;

import com.rtg.util.intervals.SequenceNameLocusSimple;

/**
 * Hold a segment comprising one or more bins in a segmentation.
 */
class Segment extends SequenceNameLocusSimple {

  private Segment mLeft;
  private Segment mRight;
  private final double mDeltaEnergy;
  private final double mSum;
  private final double mSumSq;
  private final long mBins; // Total number of original bins in this segment
  private final double mDistanceToPrevious; // Note this is the distance between midpoints
  private final int mFirstBinLength;
  private final int mLastBinLength;
  private final double mSumDistanceBetween; // mSumDistanceBetween / mBins is mean distance between bins within this segment
  private final double mSumCaseCoverage;
  private final double mSumCtrlCoverage;

  static Segment absorbRight(final Segment left, final Segment right) {
    // Discard the right segment, but adjust boundaries of the left accordingly
    if (left == null) {
      return null;
    }
    final Segment res = new Segment(left.getSequenceName(), left.getStart(), right.getEnd(),
      left.bins(), left.sum(), left.mDeltaEnergy, // Deliberately ignore energy of right
      left.distanceToPrevious(), left.mSumDistanceBetween + right.mSumDistanceBetween + right.mDistanceToPrevious,
      left.mSumCaseCoverage, left.mSumCtrlCoverage);
    res.mLeft = left.left();
    res.mRight = absorbRight(left.right(), right); // Cascade down
    return res;
  }

  static Segment absorbLeft(final Segment left, final Segment right) {
    // Discard the left segment, but adjust boundaries of the right accordingly
    if (right == null) {
      return null;
    }
    final Segment res = new Segment(left.getSequenceName(), left.getStart(), right.getEnd(),
      right.bins(), right.sum(), right.mDeltaEnergy, // Deliberately ignore energy of left
      right.distanceToPrevious(), left.mSumDistanceBetween + right.mSumDistanceBetween + right.mDistanceToPrevious,
      right.mSumCaseCoverage, right.mSumCtrlCoverage);
    res.mLeft = absorbLeft(left, right.left()); // Cascade down
    res.mRight = right.right();
    return res;
  }

  private Segment(String seqName, int start, int end, final long bins, final double sum, final double deltaEnergy,
                  final double distPrevious, final double sumDistanceBetween, final double caseCov, final double ctrlCov) {
    super(seqName, start, end);
    mFirstBinLength = end - start;
    mLastBinLength = end - start;
    mSum = sum;
    mSumSq = sum * sum;
    mBins = bins;
    mDistanceToPrevious = distPrevious;
    mSumDistanceBetween = sumDistanceBetween;
    mLeft = null;
    mRight = null;
    mDeltaEnergy = deltaEnergy;
    mSumCaseCoverage = caseCov;
    mSumCtrlCoverage = ctrlCov;
  }

  // Single bin
  Segment(String seqName, int start, int end, final double sum, final double distPrevious, double caseCov, double crtlCov) {
    this(seqName, start, end, 1, sum, Double.NEGATIVE_INFINITY, distPrevious, 0, caseCov, crtlCov);
  }

  Segment(final Segment left, final Segment right, double deltaEnergy) {
    super(left.getSequenceName(), left.getStart(), right.getEnd());
    assert left.getStart() < right.getStart() && right.getEnd() > left.getEnd();
    assert left.getSequenceName().equals(right.getSequenceName());
    mLeft = left;
    mRight = right;
    mDeltaEnergy = deltaEnergy;
    mSum = left.mSum + right.mSum;
    mSumSq = left.mSumSq + right.mSumSq;
    mBins = left.mBins + right.mBins;
    mFirstBinLength = left.mFirstBinLength;
    mLastBinLength = right.mLastBinLength;
    mDistanceToPrevious = left.mDistanceToPrevious;
    mSumDistanceBetween = left.mSumDistanceBetween + right.mSumDistanceBetween + right.mDistanceToPrevious;
    mSumCaseCoverage = left.mSumCaseCoverage + right.mSumCaseCoverage;
    mSumCtrlCoverage = left.mSumCtrlCoverage + right.mSumCtrlCoverage;
  }

  long bins() {
    return mBins;
  }

  double sum() {
    return mSum;
  }

  double sumSquares() {
    return mSumSq;
  }

  double distanceToPrevious() {
    return mDistanceToPrevious;
  }

  int firstBinLength() {
    return mFirstBinLength;
  }

  int lastBinLength() {
    return mLastBinLength;
  }

  double mean() {
    return sum() / bins();
  }

  double meanNormalizedCaseCov() {
    return mSumCaseCoverage / bins();
  }

  double meanNormalizedCtrlCov() {
    return mSumCtrlCoverage / bins();
  }

  double meanDistanceBetween() {
    return mBins == 1 ? 0 : mSumDistanceBetween / (mBins - 1);
  }

  Segment left() {
    return mLeft;
  }

  Segment right() {
    return mRight;
  }

  double deltaEnergy() {
    return mDeltaEnergy;
  }

  @Override
  public String toString() {
    return String.valueOf(bins());
  }
}
