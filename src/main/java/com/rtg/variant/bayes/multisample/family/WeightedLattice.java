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


import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.MathUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Probability weights over a lattice of subsets.
 */
@TestClass("com.rtg.variant.bayes.multisample.family.DefaultWeightedLatticeTest")
public abstract class WeightedLattice extends IntegralAbstract {

  protected final PossibilityArithmetic mArith;

  protected final BitSet mBitSet;

  /**
   * @param arith arithmetic to be used during all calculations.
   * @param bitSet methods for handling bit sets represented as integers.
   */
  WeightedLattice(PossibilityArithmetic arith, BitSet bitSet) {
    super();
    mArith = arith;
    mBitSet = bitSet;
  }

  /**
   * Used for visiting each member of the lattice which has a non-zero weight.
   */
  public interface Visitor {

    /**
     * Called once for each non-zero entry in lattice.
     * @param set for the entry.
     * @param weight for the entry.
     */
    void value(int set, double weight);
  }

  /**
   * Increment the weight associated with set.
   * @param set being incremented.
   * @param weight used for increment.
   */
  abstract void increment(int set, double weight);

  /**
   * Get the weight associated with set.
   * @param set specifies weight to be returned.
   * @return the weight.
   */
  public abstract double get(int set);

  /**
   * Assign the weight associated with set.
   * @param set being assigned.
   * @param weight used for assignment.
   */
  public abstract void set(int set, double weight);

  /**
   * Compute the lattice product.
   * @param a value being multiplied.
   * @return the product.
   */
  WeightedLattice product(final WeightedLattice a) {
    final WeightedLattice res = new DefaultWeightedLattice(mArith, mBitSet);
    final Visitor v = new Visitor() {
      @Override
      public void value(final int setv, final double valuev) {
        final Visitor w = new Visitor() {
          @Override
          public void value(final int setw, final double valuew) {
            res.increment(setv & setw, mArith.multiply(valuev, valuew)); //TODO a version for the sum to be used for recessive diseases
          }
        };
        a.visit(w);
      }
    };
    visit(v);
    return res;
  }

  /**
   * Check if two lattices are approximately equal.
   * @param that lattice being checked.
   * @param tolerance the amount of tolerance in the comparison
   * @return true iff the two lattices have all weights within tolerance of each other..
   */
  boolean approxEqual(final WeightedLattice that, final double tolerance) {
    final boolean[] equal = new boolean[1];
    equal[0] = true;
    final Visitor v = new Visitor() {
      @Override
      public void value(final int setv, final double thatv) {
        final double thisv = get(setv);
        if (!MathUtils.approxEquals(thatv, thisv, tolerance)) {
          equal[0] = false;
        }
      }
    };
    that.visit(v);
    return equal[0];
  }

  /**
   * Visit each non-zero member once.
   * @param visitor called once for each non-zero entry.
   */
  public abstract void visit(Visitor visitor);

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("[");
    final Visitor v = new Visitor() {
      int mCount = 0;
      @Override
      public void value(int set, double value) {
        if (mCount > 0) {
          sb.append(", ");
        }
        sb.append(mBitSet.toString(set)).append(":").append(Utils.realFormat(mArith.poss2Ln(value), 3));
        ++mCount;
      }
    };
    visit(v);
    sb.append("]");
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mArith);
    Exam.assertNotNull(mBitSet);
    return true;
  }

}
