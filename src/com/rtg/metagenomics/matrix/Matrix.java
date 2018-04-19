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
