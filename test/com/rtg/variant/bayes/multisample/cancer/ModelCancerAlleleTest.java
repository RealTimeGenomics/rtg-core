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

import static com.rtg.util.StringUtils.LS;

import com.rtg.reference.Ploidy;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Evidence;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.HypothesesPowerSet;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class ModelCancerAlleleTest extends TestCase {

  private static final String EXPECT_1G = "Allele Cancer Model" + LS
    + "      A  -4.248" + LS
    + "      C  -4.248" + LS
    + "    A:C  -4.248" + LS
    + "      G  -1.338" + LS
    + "    A:G  -1.984" + LS
    + "    C:G  -1.984" + LS
    + "  A:C:G  -2.345" + LS
    + "      T  -4.248" + LS
    + "    A:T  -4.248" + LS
    + "    C:T  -4.248" + LS
    + "  A:C:T  -4.248" + LS
    + "    G:T  -1.984" + LS
    + "  A:G:T  -2.345" + LS
    + "  C:G:T  -2.345" + LS
    + "A:C:G:T  -2.590" + LS;

  public void test1G() {
    final Hypotheses<Description> hyps = AbstractSomaticCallerTest.simpleHyps(0.99, 0, Ploidy.HAPLOID);
    final Description desc = hyps.description();
    final HypothesesPowerSet<Description> hypc = new HypothesesPowerSet<>(desc, SimplePossibility.SINGLETON, 0);
    final ModelCancerAllele<Description> model = new ModelCancerAllele<>(hypc, new StatisticsSnp(hypc.description()));
    model.integrity();
    final Evidence evg = new EvidenceQ(desc, 2, 0, 0, 0.05, 0.05, true, false, true, false, false);
    model.increment(evg);
    assertEquals(EXPECT_1G, model.toString());
    model.integrity();
  }
}
