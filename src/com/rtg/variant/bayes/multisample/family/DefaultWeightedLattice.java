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

import com.rtg.util.integrity.Exam;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public class DefaultWeightedLattice extends WeightedLattice {

  private final double[] mWeights;

  static WeightedLattice identity(PossibilityArithmetic arith, BitSet bitSet) {
    final DefaultWeightedLattice res = new DefaultWeightedLattice(arith, bitSet);
    res.set(bitSet.complement(0), arith.one());
    return res;
  }

  /**
   * @param arith arithmetic to be used for internal calculations.
   * @param bitSet set of bits used for subset lattice.
   */
  public DefaultWeightedLattice(PossibilityArithmetic arith, BitSet bitSet) {
    super(arith, bitSet);

    final int length = bitSet.length();
    if (length > 31) {
      throw new IllegalArgumentException("" + length);
    }
    mWeights = new double[1 << length];
    Arrays.fill(mWeights, arith.zero());
  }

  //For testing only
  DefaultWeightedLattice(final PossibilityArithmetic arith, final BitSet bitSet, final double[] values) {
    this(arith, bitSet);
    for (int i = 0; i < values.length; i++) {
      set(i, values[i]);
    }
  }

  @Override
  void increment(int set, double weight) {
    mWeights[set] = mArith.add(weight, mWeights[set]);
  }

  @Override
  public double get(int set) {
    return mWeights[set];
  }

  @Override
  public void set(int set, double weight) {
    mWeights[set] = weight;
  }

  @Override
  public void visit(Visitor visitor) {
    for (int i = 0; i < mWeights.length; i++) {
      final double w = mWeights[i];
      if (!mArith.isZero(w)) {
        visitor.value(i, w);
      }
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertEquals(1 << mBitSet.length(), mWeights.length);
    for (int i = 0; i < mWeights.length; i++) {
      Exam.assertTrue("[" + i + "]" + mWeights[i], mArith.isValidPoss(mWeights[i]));
    }
    return true;
  }

}
