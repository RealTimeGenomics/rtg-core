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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.bayes.Code;

/**
 * A non-commutative code which is a Cartesian cross product of another code.
 */
public class CodeCross implements Code {

  private final int mBaseSize;
  private final int mSize;

  /**
   * @param size of the underlying code whose product is being taken.
   */
  CodeCross(final int size) {
    assert size > 0 && size <= 46340;
    mBaseSize = size;
    mSize = mBaseSize * mBaseSize;
  }

  @Override
  public int a(int n) {
    return n / mBaseSize;
  }

  @Override
  public int b(int n) {
    return n % mBaseSize;
  }

  @Override
  public int bc(int n) {
    return b(n);
  }

  @Override
  public boolean homozygous(int n) {
    return (n / mBaseSize) == (n % mBaseSize);
  }

  @Override
  public boolean valid(int hyp) {
    return 0 <= hyp && hyp < mSize;
  }

  @Override
  public int size() {
    return mSize;
  }

  @Override
  public int rangeSize() {
    return mBaseSize;
  }

  @Override
  public int code(int a) {
    if (a >= mBaseSize) {
      throw new IllegalArgumentException(String.valueOf(a));
    }
    return a;
  }

  @Override
  public int code(int a, int b) throws IllegalArgumentException, UnsupportedOperationException {
    return a * mBaseSize + b;
  }
}
