/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.cnv.preprocess;

import java.util.Arrays;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.MathUtils;

/**
 * Holds a column of numeric data.
 */
@TestClass("com.rtg.variant.cnv.preprocess.IntColumnTest")
public class NumericColumn extends Column {

  private double[] mData = new double[100];
  private int mSize = 0;
  protected final String mFormat;

  /**
   * Construct a column with the default format
   * @param name name of the column
   */
  public NumericColumn(String name) {
    this(name, "%g");
  }

  NumericColumn(String name, String format) {
    super(name);
    mFormat = format;
  }

  @Override
  public int size() {
    return mSize;
  }

  @Override
  String toString(int i) {
    return String.format(mFormat, get(i));
  }

  @Override
  NumericColumn filter(RegionPredicate p) {
    final double[] newData = new double[mSize];
    int j = 0;
    for (int i = 0; i < mSize; ++i) {
      if (p.test(i)) {
        newData[j++] = mData[i];
      }
    }
    final NumericColumn result = new NumericColumn(getName(), mFormat);
    result.set(Arrays.copyOf(newData, j));
    return result;
  }

  @Override
  void add(String strValue) {
    double value;
    try {
      value = Double.parseDouble(strValue);
    } catch (NumberFormatException e) {
      value = Double.NaN;
    }
    add(value);
  }

  @Override
  void remove(int i) {
    if (i + 1 < mSize) {
      System.arraycopy(mData, i + 1, mData, i, 1);
    }
    mSize--;
  }

  /**
   * Append the specified value to this column.
   * @param value value to append
   */
  public void add(double value) {
    if (mSize == mData.length) {
      mData = Arrays.copyOf(mData, (int) (1 + mSize * 1.3));
    }
    mData[mSize++] = value;
  }

  void set(double[] values) {
    mSize = values.length;
    mData = Arrays.copyOf(values, values.length);
  }

  void set(StringColumn c) {
    mSize = 0;
    mData = new double[c.size()];
    for (String s : c.values()) {
      add(s);
    }
  }

  void set(int i, double value) {
    if (i >= mSize) {
      throw new ArrayIndexOutOfBoundsException();
    }
    mData[i] = value;
  }

  /**
   * Gets the value of a row
   * @param i the row
   * @return the data value
   */
  public double get(int i) {
    if (i >= mSize) {
      throw new ArrayIndexOutOfBoundsException();
    }
    return mData[i];
  }

  double[] getValues() {
    return Arrays.copyOf(mData, mSize);
  }

  double sum() {
    double sum = 0;
    for (int i = 0; i < mSize; ++i) {
      sum += mData[i];
    }
    return sum;
  }

  double mean() {
    return sum() / size();
  }

  double median() {
    return MathUtils.median(getValues());
  }
}
