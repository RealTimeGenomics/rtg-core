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

import com.rtg.mode.DNA;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Helper class for representing small sets (up to 32 items) in integers.
 */
public class BitSet extends IntegralAbstract {

  /** BitSet for nucleotides. */
  public static final BitSet DNA_SET = new BitSet(DNA.A.toString(), DNA.C.toString(), DNA.G.toString(), DNA.T.toString());

  private final int mLength;

  private final int mMask;

  private final String[] mNames;

  /**
   * @param names of values represented in sets.
   */
  BitSet(final String... names) {
    mNames = names;
    mLength = names.length;
    //be careful about case when 32 values
    final int t = 1 << (mLength - 1);
    mMask = (t - 1) + t;
  }

  /**
   * @param set of values.
   * @return human readable string with the set of values.
   */
  public String toString(final int set) {
    if (set == 0) {
      return "empty";
    }
    final StringBuilder sb = new StringBuilder();
    final Visitor vv = new Visitor() {
      int mCount = 0;
      @Override
      public void visit(final int value) {
        if (mCount > 0) {
          sb.append(":");
        }
        sb.append(mNames[value]);
        ++mCount;
      }

    };
    visit(vv, set);
    return sb.toString();
  }

  /**
   * @return the number of items that can be stored in a set.
   */
  int length() {
    return mLength;
  }

  /**
   * Construct a set given individual numbers.
   * @param values to be included in the set.
   * @return the set.
   */
  public int toSet(final int... values) {
    int t = 0;
    for (final int v : values) {
      t = t | (1 << v);
    }
    assert (t & ~mMask) == 0;
    return t;
  }

  /**
   * Visits all members of the set.
   */
  public interface Visitor {

    /**
     * Called once for each member of the set.
     * @param value the integer value of the member.
     */
    void visit(int value);
  }

  /**
   * Check if the value is contained in the set.
   * @param set to be checked.
   * @param value to be checked.
   * @return true iff value is contained in the set.
   */
  public boolean contains(final int set, final int value) {
    assert (set & ~mMask) == 0;
    assert 0 <= value && value < mLength;
    return (set & (1 << value)) != 0;
  }

  /**
   * Iterate over all values in set.
   * @param visitor call visit method on this for each value in the set.
   * @param set to be iterated over.
   */
  public void visit(final Visitor visitor, final int set) {
    int t = set;
    int i = 0;
    while (t != 0) {
      if ((t & 1) != 0) {
        visitor.visit(i);
      }
      ++i;
      t = t >>> 1;
    }
  }

  @Override
  public String toString() {
    return toString(mMask);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(0 < mLength && mLength <= Integer.SIZE);
    Exam.assertEquals(mNames.length, mLength);
    for (int i = 0; i < mLength; ++i) {
      Exam.assertFalse(mNames[i].trim().isEmpty());
      final int t = 1 << i;
      Exam.assertTrue((mMask & t) != 0);
    }
    return true;
  }

  /**
   * @param set to be complemented.
   * @return the complement of the set.
   */
  public int complement(int set) {
    return mMask & ~set;
  }

}
