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
package com.rtg.metagenomics.matrix;

import java.util.Arrays;

import com.rtg.util.Utils;

/**
 */
public class Vector {

  private final int mSize;

  private final double[] mVector;

  /**
   * For testing.
   * @param v array of values to be converted to a vector.
   */
  public Vector(final double[] v) {
    mSize = v.length;
    mVector = v.clone();
  }

  /**
   * Copy constructor.
   * @param vec the Vector to copy
   */
  public Vector(final Vector vec) {
    mSize = vec.size();
    mVector = Arrays.copyOf(vec.mVector, mSize);
  }

  /**
   * @param size length of vector.
   */
  public Vector(final int size) {
    mSize = size;
    mVector = new double[size];
  }

  /**
   * Get the <code>i'th</code> value.
   * @param i  index.
   * @return the value.
   */
  public double get(final int i) {
    return mVector[i];
  }

  /**
   * Set the <code>i'th</code> value.
   * @param i index.
   * @param v value.
   */
  public void set(final int i, final double v) {
    mVector[i] = v;
  }

  /**
   * Increment the <code>i'th</code> value.
   * @param i index.
   * @param v value.
   */
  public void incr(final int i, final double v) {
    mVector[i] += v;
  }

  /**
   * Get the length of the vector.
   * @return the length (&ge; 0).
   */
  public int size() {
    return mSize;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < mSize; ++i) {
      sb.append("  ").append(Utils.realFormat(mVector[i], 4));
    }
    return sb.toString();
  }
}
