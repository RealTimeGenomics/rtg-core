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
