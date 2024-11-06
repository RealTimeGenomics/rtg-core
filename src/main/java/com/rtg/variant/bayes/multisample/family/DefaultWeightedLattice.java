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

import com.rtg.util.integrity.Exam;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Default implementation of a weighted lattice.
 */
public final class DefaultWeightedLattice extends WeightedLattice {

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
      throw new IllegalArgumentException("length=" + length);
    }
    mWeights = new double[1 << length];
    Arrays.fill(mWeights, arith.zero());
  }

  //For testing only
  DefaultWeightedLattice(final PossibilityArithmetic arith, final BitSet bitSet, final double[] values) {
    this(arith, bitSet);
    for (int i = 0; i < values.length; ++i) {
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
    for (int i = 0; i < mWeights.length; ++i) {
      final double w = mWeights[i];
      if (!mArith.isZero(w)) {
        visitor.value(i, w);
      }
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertEquals(1 << mBitSet.length(), mWeights.length);
    for (int i = 0; i < mWeights.length; ++i) {
      Exam.assertTrue("[" + i + "]" + mWeights[i], mArith.isValidPoss(mWeights[i]));
    }
    return true;
  }

}
