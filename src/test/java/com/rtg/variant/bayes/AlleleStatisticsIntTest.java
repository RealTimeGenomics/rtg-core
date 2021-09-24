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


import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;

/**
 */
public class AlleleStatisticsIntTest extends AbstractAlleleStatisticsTest<AlleleStatisticsInt> {
  @Override
  protected AlleleStatisticsInt getAlleleStatistics(Description d) {
    return new AlleleStatisticsInt(d);
  }

  @Override
  protected void incrementAlleleStatistics(AlleleStatisticsInt stats, EvidenceInterface dist) {
    final double r = dist.mapError();
    final double q = dist.error();
    final double e = r + (1.0 - r) * q;
    stats.increment(dist, dist.read(), e);
  }

  public void testQa() {
    final AlleleStatisticsInt alleleStatistics = getAlleleStatistics(DescriptionSnp.SINGLETON);
    alleleStatistics.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0.1, 0.1, true, true, true, true, false), 1, 0.1);
    alleleStatistics.increment(new EvidenceQ(DescriptionSnp.SINGLETON, 0, 0.1, 0.1, true, true, true, true, false), 1, 0.1);
    assertEquals(20.0, alleleStatistics.qa(1));
  }
}
