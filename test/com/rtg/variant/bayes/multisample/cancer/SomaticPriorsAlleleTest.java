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
