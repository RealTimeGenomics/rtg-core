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
package com.rtg.variant.bayes.multisample.lineage;

import java.util.HashMap;
import java.util.Map;

import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.MockModel;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 * Test class.
 */
public class ModelFactorTest extends TestCase {

  public void testP() {
    final double[] p = {0.1, 0.4, 0.35, 0.15};
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hypotheses = new MockHypotheses<>(DescriptionSnp.SINGLETON, arith, true, p, 0);
    final double[] post = {0.1, 0.6, 0.25, 0.05};
    final ModelInterface<?> model = new MockModel<>(hypotheses, new StatisticsSnp(hypotheses.description()), post);
    final Variable var = new Variable("G0", hypotheses.size());
    final Factor phi = new ModelFactor(var, model);
    final Map<Variable, Integer> map = new HashMap<>();
    for (int bb = 0; bb < var.size(); ++bb) {
      map.put(var, bb);
      assertEquals(model.p(bb), phi.p(map), 1e-10);
    }
    final DefaultFactor phi2 = DefaultFactor.asDefault(phi);
    for (int bb = 0; bb < var.size(); ++bb) {
      map.put(var, bb);
      assertEquals(model.p(bb), phi2.p(map), 1e-10);
    }
  }

}
