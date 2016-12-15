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

/**
 * Energy based scorer.
 */
class EnergySegmentScorer implements SegmentScorer {
  private final double mAlpha;
  private final double mAleph;

  EnergySegmentScorer(final double alpha, final double aleph) {
    mAlpha = alpha;
    mAleph = aleph;
  }

  @Override
  public double score(final Segment a, final Segment b) {
    final double m = a.mean() - b.mean();
    final double al = a.bins();
    final double bl = b.bins();
    // SAI: Original algorithm has log(|d_1-d_2|), I added 1+ to avoid potential singularity
    return (al * bl / (al + bl)) * m * m
      + mAlpha * Math.log(1 + Math.abs(a.meanDistanceBetween() - b.meanDistanceBetween()))
      + mAleph * Math.log(1 + b.distanceToPrevious());
  }
}
