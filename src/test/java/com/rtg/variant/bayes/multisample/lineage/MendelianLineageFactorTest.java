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

import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 * Test class.
 */
public class MendelianLineageFactorTest extends TestCase {

  public void testHaploid() {
    final Variable g = new Variable("G1", 4);
    final Variable gp = new Variable("G0", 4);
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hyp = new MockHypotheses<>(new DescriptionCommon("A", "C", "G", "T"), arith, true, new double[4], 0);
    final Variable n = ForwardBackwardLineage.DE_NOVO;
    final MendelianLineageFactor m = new MendelianLineageFactor(arith, g, gp, n, 0.1, hyp, hyp);
    final Map<Variable, Integer> map = new HashMap<>();
    map.put(gp, 0);
    map.put(g, 0);
    map.put(n, 0);
    assertEquals(0.9, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(0.0, m.p(map), 1e-10);

    map.put(g, 1);
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(0.1 / 3, m.p(map), 1e-10);

    /*
    for (int parent = 0; parent < gp.size(); ++parent) {
      map.put(gp, parent);
      for (int child = 0; child < g.size(); ++child) {
        map.put(g, child);
        for (int denovo = 0; denovo < n.size(); ++denovo) {
          map.put(n, denovo);
          System.out.println(parent + " " + child + " " + denovo + " " + m.p(map));
        }
      }
    }
    */
  }

  public void testDiploid() {
    final Variable g = new Variable("G1", 10);
    final Variable gp = new Variable("G0", 10);
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hyp = new MockHypotheses<>(new DescriptionCommon("A", "C", "G", "T"), arith, false, new double[10], 0);
    final Variable n = ForwardBackwardLineage.DE_NOVO;
    final double zeta = 0.1;
    final double mu = 1 - Math.sqrt(1 - zeta);
    final MendelianLineageFactor m = new MendelianLineageFactor(arith, g, gp, n, zeta, hyp, hyp);
    final Map<Variable, Integer> map = new HashMap<>();
    final Code code = hyp.code();
    final double x = (1 - mu) / (1 + mu);
    final double y = mu / ((code.size() - 1) * (1 + mu));

    map.put(gp, code.code(0, 0));
    map.put(g, code.code(0, 0));
    map.put(n, 0);
    assertEquals(x, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(0.0, m.p(map), 1e-10);

    map.put(g, code.code(0, 1));
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(2 * y, m.p(map), 1e-10);

    map.put(g, code.code(1, 1));
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(0.0, m.p(map), 1e-10);

    map.put(g, code.code(1, 2));
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(0.0, m.p(map), 1e-10);

    map.put(gp, code.code(0, 1));
    map.put(g, code.code(0, 0));
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(y, m.p(map), 1e-10);

    map.put(g, code.code(0, 1));
    map.put(n, 0);
    assertEquals(x, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(0.0, m.p(map), 1e-10);

    map.put(g, code.code(0, 2));
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(y, m.p(map), 1e-10);

    map.put(g, code.code(1, 1));
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(y, m.p(map), 1e-10);

    map.put(g, code.code(1, 2));
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(y, m.p(map), 1e-10);

    map.put(g, code.code(2, 2));
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(0.0, m.p(map), 1e-10);

    map.put(g, code.code(2, 3));
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(0.0, m.p(map), 1e-10);

    map.put(gp, code.code(0, 2));
    map.put(g, code.code(1, 2));
    map.put(n, 0);
    assertEquals(0.0, m.p(map), 1e-10);
    map.put(n, 1);
    assertEquals(y, m.p(map), 1e-10);

    /*
    for (int parent = 0; parent < gp.size(); ++parent) {
      map.put(gp, parent);
      for (int child = 0; child < g.size(); ++child) {
        map.put(g, child);
        for (int denovo = 0; denovo < n.size(); ++denovo) {
          map.put(n, denovo);
          System.out.println(hyp.code().a(parent) + ":" + hyp.code().bc(parent) + " " + hyp.code().a(child) + ":" + hyp.code().bc(child) + " " + denovo + " " + m.p(map));
        }
      }
    }
    */
  }
  public void testAsDefault() {
    final Variable g = new Variable("G1", 10);
    final Variable gp = new Variable("G0", 10);
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final MockHypotheses<DescriptionCommon> hyp = new MockHypotheses<>(new DescriptionCommon("A", "C", "G", "T"), arith, false, new double[10], 0);
    final Variable n = ForwardBackwardLineage.DE_NOVO;
    final double zeta = 0.1;
    final MendelianLineageFactor m = new MendelianLineageFactor(arith, g, gp, n, zeta, hyp, hyp);
    final Map<Variable, Integer> map = new HashMap<>();
    final DefaultFactor def = m.asDefault();
    for (int i = 0; i < g.size(); ++i) {
      map.put(g, i);
      for (int j = 0; j < gp.size(); ++j) {
        map.put(gp, j);
        for (int k = 0; k < n.size(); ++k) {
          map.put(n, k);
          assertEquals(m.p(map), def.p(map));
        }
      }
    }
  }

}
