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
    for (int i = 0; i < dimension; ++i) {
      mMatrix[i] = new double[i + 1];
    }
  }

  @Override
  public double get(final int row, final int col) {
    if (col > row) {
      return mMatrix[col][row];
    }
    return mMatrix[row][col];
  }

  @Override
  public void set(final int row, final int col, final double v) {
    if (col > row) {
      mMatrix[col][row] = v;
    } else {
      mMatrix[row][col] = v;
    }
  }

  @Override
  public void incr(final int row, final int col, final double v) {
    if (col > row) {
      mMatrix[col][row] += v;
    } else {
      mMatrix[row][col] += v;
    }
  }

  @Override
  public int rows() {
    return mSize;
  }

  @Override
  public int cols() {
    return mSize;
  }

  @Override
  public boolean isSymmetric() {
    return true;
  }

  @Override
  public Jama.Matrix toJama() {
    final Jama.Matrix m = new Jama.Matrix(rows(), rows());
    for (int i = 0; i < rows(); ++i) {
      for (int j = 0; j < rows(); ++j) {
        m.set(i, j, get(i, j));
      }
    }
    return m;
  }
}
