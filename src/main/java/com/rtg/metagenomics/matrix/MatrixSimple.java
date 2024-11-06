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
