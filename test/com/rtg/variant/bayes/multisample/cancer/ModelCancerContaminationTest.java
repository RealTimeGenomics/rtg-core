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

package com.rtg.variant.bayes.multisample.cancer;

import static com.rtg.util.StringUtils.LS;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Evidence;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class ModelCancerContaminationTest extends TestCase {

  private HypothesesCancer<Hypotheses<Description>> getTestCats() {
    final Hypotheses<Description> hyps = AbstractSomaticCallerTest.simpleHomoHyps(0.99, 0);
    return new HypothesesCancer<>(hyps, SimplePossibility.SINGLETON);
  }

  private static final String EXPECT_1G = ""
      + "Contaminated Cancer Model contamination=   0.200" + LS
      + "A:A  -3.564  -3.564  -0.304  -3.564" + LS
      + "C:A  -3.564  -3.564  -0.304  -3.564" + LS
      + "G:A  -1.581  -1.581  -0.089  -1.581" + LS
      + "T:A  -3.564  -3.564  -0.304  -3.564" + LS
      ;
  public void test1G() {
    final HypothesesCancer<Hypotheses<Description>> hypc = getTestCats();
    final Description desc = hypc.subHypotheses().description();
    final ModelCancerContamination<Hypotheses<Description>> model = new ModelCancerContamination<>(hypc, 0.2, new StatisticsSnp(hypc.description()));
    model.integrity();
    final Evidence evg = new EvidenceQ(desc, 2, 0, 0, 0.05, 0.05, true, false, false, false);
    model.increment(evg);
    assertEquals(EXPECT_1G, model.toString());
    model.integrity();
  }
}
