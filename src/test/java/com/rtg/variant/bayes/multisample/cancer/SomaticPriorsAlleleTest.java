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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.HypothesesPowerSet;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.HypothesesCommon;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class SomaticPriorsAlleleTest extends TestCase {

  public void testHaploidNormal() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C", "G", "T");
    final HypothesesCommon<DescriptionCommon> normalHyp = new HypothesesCommon<>(desc, SimplePossibility.SINGLETON, true, 0);
    final HypothesesPowerSet<DescriptionCommon> cancerHyp = new HypothesesPowerSet<>(desc, SimplePossibility.SINGLETON, 0);
    final double[][] q = SomaticPriorsAllele.makeQ(0.1, normalHyp, cancerHyp);
    // System.out.println(Arrays.deepToString(q));
    assertEquals(normalHyp.size(), q.length);
    assertEquals(cancerHyp.size(), q[0].length);
    // 0th column should be unused
    // q[k][1 << k] is the identity
    // q[k][15] means every allele selected, should be small
    for (int k = 0; k < q.length; ++k) {
      assertEquals(0.733, q[k][(1 << k) - 1], 1e-3);
      assertEquals(7.33E-4, q[k][q[k].length - 1], 1e-7);
    }
  }

  public void testDiploidNormal() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C", "G", "T");
    final HypothesesCommon<DescriptionCommon> normalHyp = new HypothesesCommon<>(desc, SimplePossibility.SINGLETON, false, 0);
    final HypothesesPowerSet<DescriptionCommon> cancerHyp = new HypothesesPowerSet<>(desc, SimplePossibility.SINGLETON, 0);
    final double[][] q = SomaticPriorsAllele.makeQ(0.1, normalHyp, cancerHyp);
    //System.out.println(Arrays.deepToString(q));
    assertEquals(normalHyp.size(), q.length);
    assertEquals(cancerHyp.size(), q[0].length);
    // 0th column should be unused
    // q[k][(1 << ka) | (1 << kb)] is the identity, results depends on number of normal alleles
    // q[k][15] means every allele selected, should be small but result depends on number of normal alleles
    final Code code = normalHyp.code();
    for (int k = 0; k < q.length; ++k) {
      final int alleleBits = (1 << code.a(k)) | (1 << code.bc(k));
      assertEquals(Integer.bitCount(alleleBits) == 1 ? 0.733 : 0.6877, q[k][alleleBits - 1], 1e-3);
      assertEquals(Integer.bitCount(alleleBits) == 1 ? 7.33E-4 : 6.877E-3, q[k][q[k].length - 1], 1e-6);
    }
  }
}
