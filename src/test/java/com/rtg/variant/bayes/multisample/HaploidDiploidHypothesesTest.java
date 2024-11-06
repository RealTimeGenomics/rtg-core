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

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.multisample.population.DescriptionCounts;
import com.rtg.variant.bayes.snp.DescriptionNone;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class HaploidDiploidHypothesesTest extends TestCase {

  public void testInitialization() {
    final HypothesesPrior<Description> hap = new HypothesesPrior<>(DescriptionNone.SINGLETON, SimplePossibility.SINGLETON, new double[]{}, true, 0);
    final HypothesesPrior<Description> dip = new HypothesesPrior<>(DescriptionNone.SINGLETON, SimplePossibility.SINGLETON, new double[]{}, false, 1);
    final DescriptionCounts dc = new DescriptionCounts(2, 0);
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hap, dip, true, dc);
    assertEquals(hap, hdh.haploid());
    assertEquals(dip, hdh.diploid());
    assertTrue(hdh.isDefault());
    assertNotNull(hdh.getDescriptionCounts());
  }
}
