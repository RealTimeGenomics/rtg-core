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
    try {
      final NumericColumn result = (NumericColumn) clone();
      result.set(Arrays.copyOf(newData, j));
      return result;
    } catch (CloneNotSupportedException e) {
      throw new RuntimeException(e);
    }
  }

  @Override
  void add(String strValue) {
    add(toDouble(strValue));
  }

  @Override
  void add(int i, String strValue) {
    add(i, toDouble(strValue));
  }

  protected double toDouble(String strValue) {
    double value;
    try {
      value = Double.parseDouble(strValue);
    } catch (NumberFormatException e) {
      value = Double.NaN;
    }
    return value;
  }

  @Override
  void remove(int i) {
    if (i + 1 < mSize) {
      System.arraycopy(mData, i + 1, mData, i, mSize - (i + 1));
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

  /**
   * Add a value at a specific index
   * @param i the index of the entry to insert at
   * @param value the value to add
   */
  public void add(int i, double value) {
    if (i > mSize) {
      throw new ArrayIndexOutOfBoundsException();
    }
    if (mSize == mData.length) {
      mData = Arrays.copyOf(mData, (int) (1 + mSize * 1.3));
    }
    if (i < mSize) {
      System.arraycopy(mData, i, mData, i + 1, mSize - i);
    }
    mData[i] = value;
    mSize++;
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
