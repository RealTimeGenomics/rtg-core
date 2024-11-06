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

package com.rtg.variant.bayes.multisample;

import com.rtg.variant.bayes.AlleleStatisticsInt;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.EvidenceQ;

import junit.framework.TestCase;

/**
 */
public class VariantAlleleTriggerTest extends TestCase {
  static final VariantAlleleTrigger TRIGGER = new VariantAlleleTrigger(3, 0.5);
  static final DescriptionCommon DESCRIPTION = new DescriptionCommon("X", "Y", "Z");

  void evidence(AlleleStatisticsInt stats, int[] counts) {
    for (int read = 0; read < counts.length; ++read) {
      for (int count = 0; count < counts[read]; ++count) {
        stats.increment(new EvidenceQ(stats.getDescription(), read, 0.1, 0.1, true, true, true, true, true), read, 0.1);
      }
    }
  }
  public void testVariantAlleleFound() {
    final AlleleStatisticsInt stats = new AlleleStatisticsInt(DESCRIPTION);
    evidence(stats, new int[] {1, 4});
    final int variantAllele = TRIGGER.getVariantAllele(stats, DESCRIPTION, "X");
    assertEquals(1, variantAllele);

  }
  public void testBelowBothThresholds() {
    final AlleleStatisticsInt stats = new AlleleStatisticsInt(DESCRIPTION);
    evidence(stats, new int[] {2, 1});
    final int variantAllele = TRIGGER.getVariantAllele(stats, DESCRIPTION, "X");
    assertEquals(-1, variantAllele);

  }
  public void testBelowCountThreshold() {
    final AlleleStatisticsInt stats = new AlleleStatisticsInt(DESCRIPTION);
    evidence(stats, new int[] {1, 1});
    final int variantAllele = TRIGGER.getVariantAllele(stats, DESCRIPTION, "X");
    assertEquals(-1, variantAllele);
  }
  public void testBelowFractionTreshold() {
    final AlleleStatisticsInt stats = new AlleleStatisticsInt(DESCRIPTION);
    evidence(stats, new int[] {9, 4});
    final int variantAllele = TRIGGER.getVariantAllele(stats, DESCRIPTION, "X");
    assertEquals(-1, variantAllele);
  }
}
