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
