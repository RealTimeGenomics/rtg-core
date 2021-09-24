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

package com.rtg.variant.bayes.multisample.family;

import java.util.Arrays;

import com.rtg.variant.util.arithmetic.PossibilityArithmetic;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class DefaultWeightedLatticeTest extends TestCase {

  public void testIdentity() {
    final WeightedLattice identity = DefaultWeightedLattice.identity(SimplePossibility.SINGLETON, BitSet.DNA_SET);
    assertEquals("[A:C:G:T:0.000]", identity.toString());
  }

  public void test() {
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final double zero = arith.zero();
    final WeightedLattice wl = new DefaultWeightedLattice(arith, new BitSet("a", "b", "c"), new double[] {zero, arith.prob2Poss(0.2), arith.prob2Poss(0.5), arith.prob2Poss(0.4), arith.prob2Poss(0.3)});
    wl.globalIntegrity();
    assertEquals("[a:-1.609, b:-0.693, a:b:-0.916, c:-1.204]", wl.toString());
    final WeightedLattice pr = wl.product(wl);
    pr.globalIntegrity();
    assertEquals("[empty:-0.151, a:-1.609, b:-0.431, a:b:-1.833, c:-2.408]", pr.toString());
  }

  public void testApproxEqual() {
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final double zero = arith.zero();
    final WeightedLattice wla = new DefaultWeightedLattice(arith, new BitSet("a", "b", "c"), new double[] {zero, arith.prob2Poss(0.2), arith.prob2Poss(0.5), arith.prob2Poss(0.4), arith.prob2Poss(0.3)});
    wla.globalIntegrity();
    assertTrue(wla.approxEqual(wla, 0.0));
    final WeightedLattice wlb = new DefaultWeightedLattice(arith, new BitSet("a", "b", "c"), new double[] {zero, arith.prob2Poss(0.201), arith.prob2Poss(0.5), arith.prob2Poss(0.4), arith.prob2Poss(0.3)});
    wlb.globalIntegrity();
    assertTrue(wlb.approxEqual(wlb, 0.001));
    assertTrue(wla.approxEqual(wlb, 0.001));
    assertTrue(wlb.approxEqual(wla, 0.001));
  }

  public void test0() {
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final WeightedLattice wl = new DefaultWeightedLattice(arith, new BitSet("a", "b", "c"), new double[] {arith.prob2Poss(0.2)});
    wl.globalIntegrity();
    assertEquals("[empty:-1.609]", wl.toString());
  }

  public void testTooBig() {
    final PossibilityArithmetic arith = SimplePossibility.SINGLETON;
    final double[] longTest = new double[32];
    Arrays.fill(longTest, arith.zero());
    final String[] names = new String[32];
    for (int i = 0; i < 32; ++i) {
      names[i] = "" + i;
    }
    try {
      new DefaultWeightedLattice(arith, new BitSet(names), longTest);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("length=32", e.getMessage());
    }
  }
}
