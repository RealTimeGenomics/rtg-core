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
      throw new IllegalArgumentException("" + a);
    }
    return a;
  }

  @Override
  public int code(int a, int b) throws IllegalArgumentException, UnsupportedOperationException {
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
