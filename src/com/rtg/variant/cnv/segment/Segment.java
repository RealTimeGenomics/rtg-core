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

  private final double mSum;
  private final double mSumSq;
  private final long mBins; // Total number of original bins in this segment
  private final double mDistanceToPrevious; // Note this is the distance between midpoints
  private final int mFirstBinLength;
  private final int mLastBinLength;
  private final double mSumDistanceBetween; // d of VEGAWES, mSumDistanceBetween / mBins is mean distance between bins within this segment

  private Segment(int start, int end, double sum, double sumSquares, long bins, int firstBinLength, int lastBinLength, double distanceToPrevious, double distanceBetween) {
    super(start, end);
    mFirstBinLength = firstBinLength;
    mLastBinLength = lastBinLength;
    mSum = sum;
    mSumSq = sumSquares;
    mBins = bins;
    mDistanceToPrevious = distanceToPrevious;
    mSumDistanceBetween = distanceBetween;
  }

  // Single bin
  Segment(int start, int end, final double sum, final double distPrevious) {
    this(start, end, sum, sum * sum, 1, end - start, end - start, distPrevious, 0);
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

  Segment merge(final Segment other) {
    assert getStart() < other.getStart() && other.getEnd() > getEnd();
    return new Segment(getStart(),
      other.getEnd(),
      mSum + other.mSum,
      mSumSq + other.mSumSq,
      mBins + other.mBins,
      mFirstBinLength, other.mLastBinLength,
      mDistanceToPrevious,
      mSumDistanceBetween + other.mSumDistanceBetween + other.mDistanceToPrevious);
  }

  @Override
  public String toString() {
    return String.valueOf(bins());
  }
}
