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


import java.util.ArrayList;
import java.util.List;

import com.rtg.reference.Ploidy;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.HypothesesPowerSet;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.LogPossibility;

/**
 */
public class AlleleSomaticCallerTest extends AbstractSomaticCallerTest<Description> {

  @Override
  protected List<ModelInterface<Description>> getModel(Ploidy ploidy, double contamination, double same) {
    final List<ModelInterface<Description>> models = new ArrayList<>();
    for (int ref = 0; ref < 4; ++ref) {
      final Hypotheses<Description> hyps = simpleHyps(0.99, ref, ploidy);
      final HypothesesPowerSet<Description> hypc = new HypothesesPowerSet<>(hyps.description(), LogPossibility.SINGLETON, ref);
      models.add(new ModelCancerAllele<>(hypc, new StatisticsSnp(hyps.description())));
    }
    return models;
  }

  @Override
  protected AbstractSomaticCaller getSomaticCaller(final Hypotheses<Description> hypotheses, final VariantParams params, final double phi, final double psi) {
    return new AlleleSomaticCaller(new AlleleSomaticPriorsFactory<>(hypotheses), new AlleleSomaticPriorsFactory<>(hypotheses), params, phi, psi);
  }

  private static final String EXPECT_ALL_SAME = "chr1\t14\t.\tA\t.\t.\tPASS\tDP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD:SS\t0:3:0.293:0.000:100:0.00:6.51:0.00:30.330:0.00:A,3,0.293:3\t0:3:0.293:0.000:34:0.00:6.51:0.00:30.330:0.00:A,3,0.293:3:0\n";

  @Override
  protected String expectAllSame() {
    return EXPECT_ALL_SAME;
  }

  private static final String EXPECT_CANCER_EQ_NORMAL = "chr1\t14\t.\tA\tC\t.\tPASS\tDP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD:SS\t1:3:0.293:0.000:56:0.00:6.51:0.00:0.000,30.330:0.00:C,3,0.293:0,3\t1:3:0.293:0.000:34:0.00:6.51:0.00:0.000,30.330:0.00:C,3,0.293:0,3:1\n";

  @Override
  protected String getExpectCancerEqNormal() {
    return EXPECT_CANCER_EQ_NORMAL;
  }

  private static final String EXPECT_CANCER_EQ_REF = "chr1\t14\t.\tA\t.\t.\tPASS\tDP=6\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD:SS\t0:3:0.293:0.000:15:26.06:.:0.00:0.000:0.00:C,3,0.293:0\t0:3:0.293:0.000:15:0.00:6.51:0.00:30.330:0.00:A,3,0.293:3:0\n";

  @Override
  protected String getExpectCancerEqRef() {
    return EXPECT_CANCER_EQ_REF;
  }

  private static final String EXPECT_VAF = "chr1\t14\t.\tA\tG\t.\tPASS\tDP=3\tGT:VA:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:ADE:AD:SS:VAF\t1:.:0:0.000:.:11:.:.:.:0.000,0.000:.:.:0.0,0.0:0,0\t1:1:3:0.293:0.000:11:0.00:6.51:0.00:0.000,30.330:0.00:G,3,0.293:0.0,2.7:0,3:1:1.000\n";

  @Override
  protected String getExpectVaf() {
    return EXPECT_VAF;
  }

  private static final String EXPECT_NORMAL_EQ_REF = "chr1\t14\t.\tA\tC\t.\tPASS\tNCS=19.330;DP=10\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD:SSC:SS\t0:5:0.488:0.000:56:0.00:10.86:0.00:50.550,0.000:0.00:A,5,0.488:5,0\t0/1:5:0.488:0.000:14:10.86:.:0.00:0.000,50.550:0.00:C,5,0.488:0,5:1.4:2\n";

  @Override
  protected String getExpectNormalEqRef() {
    return EXPECT_NORMAL_EQ_REF;
  }

  private static final String EXPECT_ALL_DIFFERENT = "chr1\t14\t.\tA\tC,G\t.\tPASS\tNCS=19.083;DP=10\tGT:DP:RE:AR:GQ:ABP:SBP:RPB:AQ:PUR:RS:AD:SSC:SS\t1:5:0.488:0.000:32:0.00:10.86:0.00:0.000,50.550,0.000:0.00:C,5,0.488:0,5,0\t1/2:5:0.488:0.000:14:10.86:.:0.00:0.000,0.000,50.550:0.00:G,5,0.488:0,0,5:1.4:2\n";

  @Override
  protected String getExpectAllDifferent() {
    return EXPECT_ALL_DIFFERENT;
  }


}
