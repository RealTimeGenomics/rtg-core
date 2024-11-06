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
    for (int i = start; i <= mSize - left; ++i) {
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
