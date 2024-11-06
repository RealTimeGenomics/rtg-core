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
