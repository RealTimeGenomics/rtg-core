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
 * Defines a scoring function on segments.
 */
interface SegmentScorer {
  /**
   * Compute a score for the combination of two segments.
   *
   * @param a first segment
   * @param b second segment
   * @return score, bigger means the benefit of combining is bigger
   */
  double score(final Segment a, final Segment b);
}
