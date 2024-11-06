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
      .siteSpecificSomaticPriors(new ReferenceRanges<>(true))
      .contaminationBasis(42)
      .create();
    assertTrue(somaticParams.includeGainOfReference());
    assertTrue(somaticParams.includeGermlineVariants());
    assertEquals(0.987, somaticParams.lohPrior());
    assertEquals(0.345, somaticParams.somaticRate());
    assertEquals(42, somaticParams.contaminationBasis());
    assertNotNull(somaticParams.siteSpecificSomaticPriors());

  }
}
