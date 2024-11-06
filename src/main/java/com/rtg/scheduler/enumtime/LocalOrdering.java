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

package com.rtg.scheduler.enumtime;

import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Set;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Deals with local ordering of an overall ordering on an integer time and an enumeration.
 * @param <E> enumeration type.
 */
public class LocalOrdering<E extends Enum<E>> extends IntegralAbstract {

  //Floyd Warshall
  static void transitiveClosure(boolean[][] matrix) {
    final int n = matrix.length;
    for (int k = 0; k < n; ++k) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          matrix[i][j] |= matrix[i][k] & matrix[k][j];
        }
      }
    }
  }

  private final Set<EnumTimeId<E>>[] mTo;
  private final Set<EnumTimeId<E>>[] mFrom;
  private final boolean[][] mBefore;
  private boolean mFrozen = false;

  @SuppressWarnings("unchecked")
  private Set<EnumTimeId<E>>[] getArray(int len) {
    return (Set<EnumTimeId<E>>[]) new Set<?>[len];
  }

  /**
   * @param values array of all values in the enumeration.
   */
  public LocalOrdering(final E[] values) {
    final int length = values.length;
    mTo = getArray(length);
    mFrom = getArray(length);
    mBefore = new boolean[length][length];

    for (int i = 0; i < mTo.length; ++i) {
      mTo[i] = new HashSet<>();
    }
    for (int i = 0; i < mFrom.length; ++i) {
      mFrom[i] = new LinkedHashSet<>();
    }
  }

  /**
   * Given the index (an enumeration) return the set of time, enumeration pairs that
   * it will contribute its result to.
   * @param index the enumeration index selecting the id at a particular time.
   * @return the set of destinations for results.
   */
  public Set<EnumTimeId<E>> to(final E index) {
    if (!mFrozen) {
      throw new RuntimeException();
    }
    return mTo[index.ordinal()];
  }

  /**
   * Given the index (an enumeration) return the set of time, enumeration pairs that
   * it will draw its argument from.
   * @param index the enumeration index selecting the id at a particular time.
   * @return the set of arguments (in the predefined order that the corresponding <code>Job</code> will be expecting its arguments).
   */
  public Set<EnumTimeId<E>> from(final E index) {
    if (!mFrozen) {
      throw new RuntimeException();
    }
    return mFrom[index.ordinal()];
  }

  /**
   * Test if a is before b in the irreflexive pre-order.
   * @param a first argument.
   * @param b second argument.
   * @return true iff a is before b in the pre-order.
   */
  public boolean before(final E a, final E b) {
    if (!mFrozen) {
      throw new RuntimeException();
    }
    return mBefore[a.ordinal()][b.ordinal()];
  }

  /**
   * Freeze this object so that no more <code>setLink</code> calls can be done.
   */
  public void freeze() {
    assert !mFrozen;
    transitiveClosure(mBefore);
    mFrozen = true;
  }

  boolean frozen() {
    return mFrozen;
  }

  /**
   * Set the dependency links both for from and to sets.
   * @param from source for any transfer of results.
   * @param increment change in chunk level (usually 1 or 0).
   * @param to destination of any transfer of results.
   */
  public void setLink(final E from, final int increment, final E to) {
    if (mFrozen) {
      throw new RuntimeException();
    }
    mTo[from.ordinal()].add(new EnumTimeId<>(increment, to));
    mFrom[to.ordinal()].add(new EnumTimeId<>(-increment, from));
    if (increment == 0) {
      mBefore[from.ordinal()][to.ordinal()] = true;
    }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    if (mFrozen) {
      for (int i = 0; i < mBefore.length; ++i) {
        Exam.assertFalse(mBefore[i][i]);
      }
    }
    for (boolean[] aMBefore : mBefore) {
      Exam.assertEquals(mBefore.length, aMBefore.length);
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertEquals(mBefore.length, mTo.length);
    Exam.assertEquals(mBefore.length, mFrom.length);
    return true;
  }

}
