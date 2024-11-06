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

import static com.rtg.util.StringUtils.LS;

import com.rtg.reference.Ploidy;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Evidence;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class ModelCancerContaminationTest extends TestCase {

  private HypothesesCancer<Hypotheses<Description>> getTestCats() {
    final Hypotheses<Description> hyps = AbstractSomaticCallerTest.simpleHyps(0.99, 0, Ploidy.HAPLOID);
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
    final ModelCancerContamination<Hypotheses<Description>> model = new ModelCancerContamination<>(hypc, 0.2, new StatisticsSnp(hypc.description()), new NoAlleleBalance());
    model.integrity();
    final Evidence evg = new EvidenceQ(desc, 2, 0, 0, 0.05, 0.05, true, false, true, false, false);
    model.increment(evg);
    assertEquals(EXPECT_1G, model.toString());
    model.integrity();
  }
}
