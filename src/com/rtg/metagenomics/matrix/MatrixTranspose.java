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
public class MatrixTranspose extends Matrix {

  private final Matrix mMatrix;

  /**
   * @param matrix to be transposed.
   */
  public MatrixTranspose(final Matrix matrix) {
    mMatrix = matrix;
  }

  @Override
  public double get(final int row, final int col) {
    return mMatrix.get(col, row);
  }

  @Override
  public void set(final int row, final int col, final double v) {
    mMatrix.set(col, row, v);
  }

  @Override
  public void incr(final int row, final int col, final double v) {
    mMatrix.incr(col, row, v);
  }

  @Override
  public int rows() {
    return mMatrix.cols();
  }

  @Override
  public int cols() {
    return mMatrix.rows();
  }

  @Override
  public boolean isSymmetric() {
    return mMatrix.isSymmetric();
  }

  @Override
  public Jama.Matrix toJama() {
    final Jama.Matrix m = new Jama.Matrix(mMatrix.rows(), mMatrix.rows());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < rows(); ++j) {
        m.set(i, j, get(i, j));
      }
    }
    return m;
  }
}
