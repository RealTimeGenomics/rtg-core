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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.Utils;

/**
 * Two dimensional square real valued matrices.
 */
@TestClass(value = {"com.rtg.metagenomics.matrix.MatrixSimpleTest"})
public abstract class Matrix {

  /**
   * Get the <code>(row, col)th</code> value.
   * @param row first index.
   * @param col second index.
   * @return the value.
   */
  public abstract double get(int row, int col);

  /**
   * Set the <code>(row, col)th</code> value.
   * @param row first index.
   * @param col second index.
   * @param v value.
   */
  public abstract void set(int row, int col, double v);

  /**
   * Increment the <code>(row, col)th</code> value.
   * @param row first index.
   * @param col second index.
   * @param v value.
   */
  public abstract void incr(int row, int col, double v);

  /**
   * Number of rows in the matrix.
   * @return the length (&ge; 0).
   */
  public abstract int rows();

  /**
   * Number of rows in the matrix.
   * @return the length (&ge; 0).
   */
  public abstract int cols();

  /**
   * Test if is symmetric.
   * @return true iff symmetric.
   */
  public abstract boolean isSymmetric();

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < rows(); ++i) {
      sb.append("[").append(i).append("]");
      final int hi = isSymmetric() ? i : rows() - 1;
      for (int j = 0; j <= hi; ++j) {
        sb.append("  ").append(Utils.realFormat(get(i, j), 4));
      }
      sb.append(System.lineSeparator());
    }
    return sb.toString();
  }

  /**
   * Convert to <code>Jama matrix</code>
   * @return the matrix
   */
  public abstract Jama.Matrix toJama();

}
