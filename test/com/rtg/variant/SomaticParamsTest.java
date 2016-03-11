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

package com.rtg.variant;

import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class SomaticParamsTest extends TestCase {
  public void testOmnes() {
    new TestParams(SomaticParams.class, SomaticParamsBuilder.class).check();
  }

  public void test() {
    final SomaticParams somaticParams = new SomaticParamsBuilder()
      .includeGainOfReference(true)
      .includeGermlineVariants(true)
      .lohPrior(0.987)
      .somaticRate(0.345)
      .siteSpecificSomaticPriors(new ReferenceRanges<Double>(true))
      .create();
    assertTrue(somaticParams.includeGainOfReference());
    assertTrue(somaticParams.includeGermlineVariants());
    assertEquals(0.987, somaticParams.lohPrior());
    assertEquals(0.345, somaticParams.somaticRate());
    assertNotNull(somaticParams.siteSpecificSomaticPriors());

  }
}
