/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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
 * Power set without empty set.
 */
public class CodePowerSet implements Code {

  private final int mSize;

  /**
   * Code representing the power set of specified number of items but excluding empty set.
   * @param n underlying description size
   */
  public CodePowerSet(final int n) {
    assert n <= 31; // Actually n > 6 or so will make for slow calling ...
    mSize = (1 << n) - 1;
  }

  @Override
  public int a(final int n) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int b(final int n) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int bc(final int n) {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean homozygous(final int n) {
    return ((n + 1) & n) == 0;
  }

  @Override
  public boolean valid(final int hyp) {
    return 0 <= hyp && hyp < mSize;
  }

  @Override
  public int size() {
    return mSize;
  }

  @Override
  public int rangeSize() {
    throw new UnsupportedOperationException();
  }

  @Override
  public int code(final int a) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int code(final int a, final int b) {
    throw new UnsupportedOperationException();
  }
}
