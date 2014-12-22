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
package com.rtg.util.array.zeroindex;

import com.rtg.util.array.AbstractIndex;
import com.rtg.util.array.ExtensibleIndex;
import com.rtg.util.format.FormatInteger;

/**
 * All positions in the index return the same value. That is this uses zero bits per entry.
 *
 */
public final class ZeroIndex extends AbstractIndex implements ExtensibleIndex {

  private static final FormatInteger FORMAT_INTEGER = new FormatInteger(5);

  private final long mConstant;

  /**
   * @param length number of items.
   * @param constant value stored at each location.
   * @exception NegativeArraySizeException if length less than 0
   */
  public ZeroIndex(final long length, final long constant) {
    super(length);
    mConstant = constant;
    assert globalIntegrity();
  }

  @Override
  public long bytes() {
    return 0L;
  }

  @Override
  public long get(final long index) {
    check(index);
    return mConstant;
  }

  @Override
  public void set(final long index, final long value) {
    assert value == mConstant; // : value + ":" + mConstant;
  }

  @Override
  protected FormatInteger formatValue() {
    return FORMAT_INTEGER;
  }

  @Override
  public long extendBy(long increment) {
    assert increment >= 0;
    final long res = mLength;
    mLength += increment;
    return res;
  }

  @Override
  public void extendTo(long size) {
    mLength = Math.max(size, mLength);
  }

  @Override
  public void trim(long length) {
    assert length >= 0;
    mLength = length;
  }

  @Override
  public long getSigned(long offset) {
    check(offset);
    return mConstant;
  }

  @Override
  public void setSigned(long offset, long value) {
    assert value == mConstant; // : value + ":" + mConstant;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("Index [").append(length()).append("]").append(LS);
    if (mConstant != 0) {
      sb.append(mConstant).append(" constant").append(LS);
    }
  }

  @Override
  public boolean safeFromWordTearing() {
    return true;
  }
}
