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

import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

import junit.framework.TestCase;

/**
 */
public class SFunctionTest extends TestCase {

  public void test3children() {
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final BitSet bitSet = BitSet.DNA_SET;
    final WeightedLattice init = new DefaultWeightedLattice(arith, bitSet);
    final double one = arith.one();
    init.set(3, one);
    final WeightedLattice[] children = new DefaultWeightedLattice[3];
    children[0] = new DefaultWeightedLattice(arith, bitSet);
    children[0].set(3, one);
    children[1] = new DefaultWeightedLattice(arith, bitSet);
    children[1].set(6, one);
    children[2] = new DefaultWeightedLattice(arith, bitSet);
    children[2].set(12, one);
    final SFunction sf = new SFunction(init, children);
    assertEquals("[empty:0.000]", sf.all().toString());
    assertEquals("[empty:0.000]", sf.excludeChild(0).toString());
    assertEquals("[empty:0.000]", sf.excludeChild(1).toString());
    assertEquals("[C:0.000]", sf.excludeChild(2).toString());
  }

  public void test2children() {
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final BitSet bitSet = BitSet.DNA_SET;
    final WeightedLattice init = new DefaultWeightedLattice(arith, bitSet);
    final double one = arith.one();
    init.set(3, one);
    final WeightedLattice[] children = new DefaultWeightedLattice[2];
    children[0] = new DefaultWeightedLattice(arith, bitSet);
    children[0].set(3, one);
    children[1] = new DefaultWeightedLattice(arith, bitSet);
    children[1].set(6, one);
    final SFunction sf = new SFunction(init, children);
    assertEquals("[C:0.000]", sf.all().toString());
    assertEquals("[C:0.000]", sf.excludeChild(0).toString());
    assertEquals("[A:C:0.000]", sf.excludeChild(1).toString());
  }

  public void test1children() {
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final BitSet bitSet = BitSet.DNA_SET;
    final WeightedLattice init = new DefaultWeightedLattice(arith, bitSet);
    final double one = arith.one();
    init.set(3, one);
    final WeightedLattice[] children = new DefaultWeightedLattice[1];
    children[0] = new DefaultWeightedLattice(arith, bitSet);
    children[0].set(1, one);
    final SFunction sf = new SFunction(init, children);
    assertEquals("[A:0.000]", sf.all().toString());
    assertEquals("[A:C:0.000]", sf.excludeChild(0).toString());
  }

  public void test0children() {
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final BitSet bitSet = BitSet.DNA_SET;
    final WeightedLattice init = new DefaultWeightedLattice(arith, bitSet);
    final double one = arith.one();
    init.set(3, one);
    final WeightedLattice[] children = new DefaultWeightedLattice[0];
    final SFunction sf = new SFunction(init, children);
    assertEquals("[A:C:0.000]", sf.all().toString());
  }

}
