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
package com.rtg.variant.bayes;

/**
 */
public class AlleleStatisticsDoubleTest extends AbstractAlleleStatisticsTest<AlleleStatisticsDouble> {

  @Override
  protected AlleleStatisticsDouble getAlleleStatistics(Description d) {
    return new AlleleStatisticsDouble(d);
  }

  @Override
  protected void incrementAlleleStatistics(AlleleStatisticsDouble stats, EvidenceInterface dist) {
    final double r = dist.mapError();
    final double q = dist.error();
    final double e = r + (1.0 - r) * q;
    stats.increment(dist, dist.read(), e, 1.0);
  }
}
