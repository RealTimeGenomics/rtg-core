/*
 * Copyright (c) 2016. Real Time Genomics Limited.
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

/**
 * Hold a pair of adjacent segments and thier score.
 */
class AdjacentSegments {
  private final double mScore;
  private final Segment mFirst;
  private final Segment mSecond;

  /**
   * @param score the score assigned to this pairing of segments
   * @param first the first of two segments. Usually the one with lowest starting position.
   * @param second the second of two segments. Usually the higher starting position.
   */
  AdjacentSegments(double score, Segment first, Segment second) {
    this.mScore = score;
    this.mFirst = first;
    this.mSecond = second;
  }

  /** @return the score assigned to these segments */
  public double getScore() {
    return mScore;
  }

  /** @return the first of the two segments */
  public Segment getFirst() {
    return mFirst;
  }

  /** @return the second of the two segments */
  public Segment getSecond() {
    return mSecond;
  }
}
