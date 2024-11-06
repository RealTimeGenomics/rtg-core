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
package com.rtg.metagenomics.metasnp;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

import junit.framework.TestCase;

/**
 */
public class EmIterateTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void checkSimple(final boolean updateBeta) {
    final List<Integer> ref = Arrays.asList(0, 1, 2, 0);

    final List<double[][]> evidence = new ArrayList<>();
    evidence.add(new double[][] {{0, 0, 100, 200}, {0, 0, 200, 100}});
    evidence.add(new double[][] {{100, 0, 0, 200}, {200, 0, 0, 100}});
    evidence.add(new double[][] {{300, 0, 0, 0}, {300, 0, 0, 0}});
    evidence.add(new double[][] {{0, 100, 0, 200}, {0, 200, 0, 100}});
    final int strains = 2;
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final List<EmIterate.EmResult> iterations = EmIterate.iterate(ref, evidence, strains, updateBeta ? EmIterate.BetaType.REESTIMATE : EmIterate.BetaType.STATIC, 0.001, MetaSnpCli.initXi(2, 2, arith));
    final EmIterate.EmResult result = iterations.get(iterations.size() - 1);
    final double[][] xi = result.mXi;
    assertEquals(0.33, arith.poss2Prob(xi[0][0]), 0.01);
    assertEquals(0.66, arith.poss2Prob(xi[0][1]), 0.01);
    assertEquals(0.66, arith.poss2Prob(xi[1][0]), 0.01);
    assertEquals(0.33, arith.poss2Prob(xi[1][1]), 0.01);
    final int[][] expected = {{2, 3}, {0, 3}, {0, 0}, {1, 3}};
    for (int i = 0; i < result.mAssignments.size(); ++i) {
      assertTrue(Arrays.equals(expected[i], result.mAssignments.get(i).mCalls));
    }
  }
  
  public void testSimpleNoBeta() {
    checkSimple(false);
  }
  public void testSimpleBeta() {
    checkSimple(true);
  }
  
  public void testNoisy() {
    final List<Integer> ref = Arrays.asList(0, 1, 2, 0);

    final List<double[][]> evidence = new ArrayList<>();
    evidence.add(new double[][] {{2, 0, 100, 200}, {0, 1, 200, 100}});
    evidence.add(new double[][] {{100, 3, 0, 200}, {200, 0, 2, 100}});
    evidence.add(new double[][] {{300, 0, 1, 0}, {300, 0, 2, 0}});
    evidence.add(new double[][] {{0, 100, 5, 200}, {2, 200, 2, 100}});
    final int strains = 2;
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final List<EmIterate.EmResult> iterations = EmIterate.iterate(ref, evidence, strains, EmIterate.BetaType.STATIC, 0.001, MetaSnpCli.initXi(2, 2, arith));
    final EmIterate.EmResult result = iterations.get(iterations.size() - 1);
    final double[][] xi = result.mXi;
    assertEquals(0.33, arith.poss2Prob(xi[0][0]), 0.01);
    assertEquals(0.66, arith.poss2Prob(xi[0][1]), 0.01);
    assertEquals(0.66, arith.poss2Prob(xi[1][0]), 0.01);
    assertEquals(0.33, arith.poss2Prob(xi[1][1]), 0.01);
    final int[][] expected = {{2, 3}, {0, 3}, {0, 0}, {1, 3}};
    for (int i = 0; i < result.mAssignments.size(); ++i) {
      assertTrue(Arrays.equals(expected[i], result.mAssignments.get(i).mCalls));
    }

  }
}
