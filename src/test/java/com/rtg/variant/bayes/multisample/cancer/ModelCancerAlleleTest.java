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
