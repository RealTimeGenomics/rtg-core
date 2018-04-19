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
 * Real valued matrix.
 */
public class MatrixSimple extends Matrix {

  private final int mRows;
  private final int mCols;

  private final double[][] mMatrix;

  /**
   * Construct a zero matrix.
   * @param rows number of rows
   * @param cols number of columns
   */
  public MatrixSimple(final int rows, final int cols) {
    mRows = rows;
    mCols = cols;
    mMatrix = new double[rows][cols];
  }

  /**
   * Construct a zero square matrix.
   * @param dimension n of an n x n matrix.
   */
  public MatrixSimple(final int dimension) {
    this(dimension, dimension);
  }

  /**
   * Construct a simple matrix from the given matrix.
   * @param m matrix
   */
  public MatrixSimple(final Matrix m) {
    mRows = m.rows();
    mCols = m.cols();
    mMatrix = new double[mRows][mCols];
    for (int y = 0; y < mRows; ++y) {
      for (int x = 0; x < mCols; ++x) {
        mMatrix[y][x] = m.get(y, x);
      }
    }
  }

  @Override
  public double get(final int row, final int col) {
    return mMatrix[row][col];
  }

  @Override
  public void set(final int row, final int col, final double v) {
    mMatrix[row][col] = v;
  }

  @Override
  public void incr(final int row, final int col, final double v) {
    mMatrix[row][col] += v;
  }

  @Override
  public int rows() {
    return mRows;
  }

  @Override
  public int cols() {
    return mCols;
  }

  @Override
  public boolean isSymmetric() {
    if (mRows != mCols) {
      return false;
    }
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < i; ++j) {
        if (get(i, j) != get(j, i)) {
          return false;
        }
      }
    }
    return true;
  }

  @Override
  public Jama.Matrix toJama() {
    return new Jama.Matrix(mMatrix);
  }
}
