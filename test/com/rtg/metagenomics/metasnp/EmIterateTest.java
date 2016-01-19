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
    for (int i = 0; i < result.mAssignments.size(); i++) {
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
    for (int i = 0; i < result.mAssignments.size(); i++) {
      assertTrue(Arrays.equals(expected[i], result.mAssignments.get(i).mCalls));
    }

  }
}
