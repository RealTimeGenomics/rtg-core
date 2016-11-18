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

/**
 */
public class MatrixSymmetric extends Matrix {

  private final int mSize;

  private final double[][] mMatrix;

  /**
   * @param dimension n of an n x n matrix.
   */
  public MatrixSymmetric(final int dimension) {
    mSize = dimension;
    mMatrix = new double[dimension][];
    for (int i = 0; i < dimension; i++) {
      mMatrix[i] = new double[i + 1];
    }
  }

  @Override
  public double get(final int i, final int j) {
    if (j > i) {
      return mMatrix[j][i];
    }
    return mMatrix[i][j];
  }

  @Override
  public void set(final int i, final int j, final double v) {
    if (j > i) {
      mMatrix[j][i] = v;
    } else {
      mMatrix[i][j] = v;
    }
  }

  @Override
  public void incr(final int i, final int j, final double v) {
    if (j > i) {
      mMatrix[j][i] += v;
    } else {
      mMatrix[i][j] += v;
    }
  }

  @Override
  public int size() {
    return mSize;
  }

  @Override
  public boolean isSymmetric() {
    return true;
  }

  @Override
  public Jama.Matrix toJama() {
    final Jama.Matrix m = new Jama.Matrix(size(), size());
    for (int i = 0; i < size(); i++) {
      for (int j = 0; j < size(); j++) {
        m.set(i, j, get(i, j));
      }
    }
    return m;
  }
}
