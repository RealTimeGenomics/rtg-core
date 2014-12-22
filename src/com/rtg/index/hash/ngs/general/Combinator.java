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
package com.rtg.index.hash.ngs.general;

import java.util.BitSet;

/**
 */
public abstract class Combinator {

  private final int mSize;

  private final int mToBeSet;

  private final BitSet mSet;
  /**
   * @param size total number of positions.
   * @param toBeSet the number of bits to be set on.
   */
  Combinator(final int size, final int toBeSet) {
    assert 0 <= toBeSet && toBeSet <= size;
    mSize = size;
    mToBeSet = toBeSet;
    mSet = new BitSet(mSize);
  }

  /**
   * Generate all the permutations.
   */
  public void combine() {
    comb(0, mToBeSet);
  }

  void comb(final int start, final int left) {
    assert left >= 0;
    assert start <= mSize;
    if (left == 0) {
      permutation();
      return;
    }
    assert start < mSize;
    for (int i = start; i <= mSize - left; i++) {
      mSet.set(i);
      comb(i + 1, left - 1);
      mSet.clear(i);
    }
  }

  /**
   * Called once for each permutation.
   * The current state of the permutation is available via the calls to <code>size(), tobeSet() and getBit(int)</code>.
   */
  public abstract void permutation();

  protected int size() {
    return mSize;
  }

  protected int toBeSet() {
    return mToBeSet;
  }

  protected boolean getBit(final int index) {
    return mSet.get(index);
  }
}
