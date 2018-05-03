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

import com.rtg.util.intervals.Range;

/**
 * Hold a segment comprising one or more bins in a segmentation.
 */
class Segment extends Range {

  private final Segment mLeft;
  private final Segment mRight;
  private final double mSum;
  private final double mSumSq;
  private final long mBins; // Total number of original bins in this segment
  private final double mDistanceToPrevious; // Note this is the distance between midpoints
  private final int mFirstBinLength;
  private final int mLastBinLength;
  private final double mSumDistanceBetween; // mSumDistanceBetween / mBins is mean distance between bins within this segment

  // Single bin
  Segment(int start, int end, final double sum, final double distPrevious) {
    super(start, end);
    mFirstBinLength = end - start;
    mLastBinLength = end - start;
    mSum = sum;
    mSumSq = sum * sum;
    mBins = 1;
    mDistanceToPrevious = distPrevious;
    mSumDistanceBetween = 0;
    mLeft = null;
    mRight = null;
  }

  Segment(final Segment left, final Segment right) {
    super(left.getStart(), right.getEnd());
    assert left.getStart() < right.getStart() && right.getEnd() > left.getEnd();
    mLeft = left;
    mRight = right;
    mSum = left.mSum + right.mSum;
    mSumSq = left.mSumSq + right.mSumSq;
    mBins = left.mBins + right.mBins;
    mFirstBinLength = left.mFirstBinLength;
    mLastBinLength = right.mLastBinLength;
    mDistanceToPrevious = left.mDistanceToPrevious;
    mSumDistanceBetween = left.mSumDistanceBetween + right.mSumDistanceBetween + right.mDistanceToPrevious;
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

  double meanDistanceBetween() {
    return mBins == 1 ? 0 : mSumDistanceBetween / (mBins - 1);
  }

  Segment left() {
    return mLeft;
  }

  Segment right() {
    return mRight;
  }

  @Override
  public String toString() {
    return String.valueOf(bins());
  }
}
