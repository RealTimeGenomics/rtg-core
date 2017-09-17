/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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

import com.rtg.reference.Ploidy;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.util.arithmetic.LogPossibility;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class HypothesesPowerSetTest extends TestCase {

  public void test() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C", "G", "T");
    final HypothesesPowerSet<DescriptionCommon> hyp = new HypothesesPowerSet<>(desc, LogPossibility.SINGLETON, 2);
    assertEquals(15, hyp.size());
    assertFalse(hyp.valid(0));
    assertTrue(hyp.valid(1));
    assertTrue(hyp.valid(15));
    assertFalse(hyp.valid(16));
    assertEquals("A:C:G:T", hyp.name(15));
    assertEquals("G", hyp.name(hyp.reference()));
    assertEquals(LogPossibility.SINGLETON, hyp.arithmetic());
    assertEquals(7, hyp.maxNameLength());
    assertEquals(Ploidy.POLYPLOID, hyp.ploidy());
    assertTrue(hyp.code() instanceof CodePowerSet);
  }
}
