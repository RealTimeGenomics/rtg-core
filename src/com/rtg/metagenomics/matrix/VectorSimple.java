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

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class VectorSimple extends IntegralAbstract implements Vector {

  private final int mSize;

  private final double[] mVector;

  /**
   * For testing.
   * @param v array of values to be converted to a vector.
   */
  public VectorSimple(final double[] v) {
    mSize = v.length;
    mVector = v.clone();
  }

  /**
   * Copy constructor
   * @param copy the Vector to copy
   */
  public VectorSimple(Vector copy) {
    mSize = copy.dimension();
    mVector = new double[mSize];
    for (int i = 0; i < mSize; i++) {
      mVector[i] = copy.get(i);
    }
  }

  /**
   * @param size length of vector.
   */
  public VectorSimple(final int size) {
    mSize = size;
    mVector = new double[size];
  }

  @Override
  public double get(final int i) {
    return mVector[i];
  }

  @Override
  public void set(final int i, final double v) {
    mVector[i] = v;
  }

  @Override
  public void incr(final int i, final double v) {
    mVector[i] += v;
  }

  @Override
  public int dimension() {
    return mSize;
  }

  @Override
  public void toString(final StringBuilder sb) {
    for (int i = 0; i < mSize; i++) {
      sb.append("  ").append(Utils.realFormat(mVector[i], 4));
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mSize >= 0 && mVector.length == mSize);
    return true;
  }

}
