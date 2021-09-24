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
