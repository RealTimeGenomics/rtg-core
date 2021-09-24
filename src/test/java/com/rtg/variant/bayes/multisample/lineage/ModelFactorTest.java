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
