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

package com.rtg.variant.bayes;

/**
 * Map pairs of integer values a, b from the range <code>0 ... l</code> (exclusive) where a &le; b into a single
 * contiguous set of codes from <code>0 ... l * (l + 1) / 2</code> (exclusive).
 * <br><br>
 * This is done in two steps. The diagram below shows the numbering of a mathematically simpler triangular code.
 * This uses pairs <code>i, j where i + j &lt; l</code>. We use the variable <code>k = i + j</code>.
 * <Br>
 * <img src="doc-files/diagrams/trimmed1.jpg" alt="image">
 * <br><br>
 * The diagram below shows the numbering for the final a, b code which is obtained by transformation from the
 * triangular code.
 * <br>
 * <img src="doc-files/diagrams/trimmed2.jpg" alt="image">
 */
public final class CodeDiploid implements Code {
  /** Table of <code>k(n)</code> values. */
  private final int[] mValues;
  /** Table of <code>i * (i + 1) / 2</code>*/
  private final int[] mSqr;
  private final int mL;
  private final int mLs;

  /**
   * @param l number of haploid hypotheses.
   */
  public CodeDiploid(final int l) {
    mL = l - 1;
    final int len = (l * (l + 1)) / 2;
    mLs = len - 1;
    mValues = new int[len];
    mSqr = new int[l];
    int n = 0;
    for (int i = 0; i < l; ++i) {
      mSqr[i] = n;
      for (int j = 0; j <= i; ++j) {
        mValues[n] = i;
        ++n;
      }
    }
  }

  int k(final int n) {
    return mValues[n];
  }

  int s(final int k) {
    return mSqr[k];
  }

  /**
   * Extract the first part of the triangle code from n.
   * @param n to be decoded.
   * @return first part of triangle code.
   */
  int i(final int n) {
    final int k = k(n);
    final int np = s(k);
    final int j = n - np;
    return k - j;
  }

  /**
   * Extract the second part of the triangle code from n.
   * @param n to be decoded.
   * @return second part of triangle code.
   */
  int j(final int n) {
    final int k = k(n);
    final int np = s(k);
    return n - np;
  }

  @Override
  public int size() {
    return mValues.length;
  }

  @Override
  public int rangeSize() {
    return mL + 1;
  }

  @Override
  public int a(final int n) {
    final int m = mLs - n;
    final int j = j(m);
    return mL - j;
  }

  @Override
  public int b(final int n) {
    final int m = mLs - n;
    return i(m);
  }

  @Override
  public int bc(int n) {
    return b(n);
  }

  @Override
  public boolean homozygous(int n) {
    return n <= mL;
  }

  @Override
  public boolean valid(int hyp) {
    return 0 <= hyp && hyp < mValues.length;
  }

  @Override
  public int code(int a) {
    if (!valid(a)) {
      throw new IllegalArgumentException("a=" + a);
    }
    return a;
  }

  @Override
  public int code(int a, int b) {
    assert 0 <= a && a <= mL;
    if (b > a) {
      return code(b, a);
    }
    final int j = mL - a;
    final int k = b + j;
    final int x = ((k * (k + 1)) >> 1) + j;
    return mLs - x;
  }
}
