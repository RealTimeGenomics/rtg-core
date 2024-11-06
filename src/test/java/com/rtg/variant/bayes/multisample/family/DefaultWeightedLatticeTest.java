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
