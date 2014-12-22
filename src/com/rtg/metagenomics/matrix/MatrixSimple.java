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

import com.rtg.util.integrity.Exam;

/**
 */
public class MatrixSimple extends Matrix {

  private final int mSize;

  private final double[][] mMatrix;

  /**
   * Create a matrix from arrays of doubles.
   * Intended for testing.
   * Copies the array.
   * @param values to be used.
   */
  MatrixSimple(final double[][] values) {
    mSize = values.length;
    mMatrix = new double[mSize][];
    for (int i = 0; i < mSize; i++) {
      mMatrix[i] = new double[mSize];
      System.arraycopy(values[i], 0, mMatrix[i], 0, mSize);
    }
  }

  /**
   * @param dimension n of an n x n matrix.
   */
  public MatrixSimple(final int dimension) {
    mSize = dimension;
    mMatrix = new double[dimension][];
    for (int i = 0; i < dimension; i++) {
      mMatrix[i] = new double[dimension];
    }
  }

  @Override
  public double get(final int i, final int j) {
    return mMatrix[i][j];
  }

  @Override
  public void set(final int i, final int j, final double v) {
    mMatrix[i][j] = v;
  }

  @Override
  public void incr(final int i, final int j, final double v) {
    mMatrix[i][j] += v;
  }

  @Override
  public int dimension() {
    return mSize;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (int i = 0; i < mSize; i++) {
      Exam.assertEquals(mSize, mMatrix[i].length);
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mSize >= 0);
    Exam.assertFalse(isSymmetric());
    return true;
  }

  @Override
  public boolean isSymmetric() {
    return false;
  }

  @Override
  public Jama.Matrix toJama() {
    return new Jama.Matrix(mMatrix);
  }
}
