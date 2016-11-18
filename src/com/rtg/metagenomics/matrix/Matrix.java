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
   * Get the (i, j)th value.
   * @param i first index.
   * @param j second index.
   * @return the value.
   */
  public abstract double get(int i, int j);

  /**
   * Set the (i, j)th value.
   * @param i first index.
   * @param j second index.
   * @param v value.
   */
  public abstract void set(int i, int j, double v);

  /**
   * Increment the (i, j)th value.
   * @param i first index.
   * @param j second index.
   * @param v value.
   */
  public abstract void incr(int i, int j, double v);

  /**
   * Get the length n of a square n x n matrix.
   * @return the length (&ge; 0).
   */
  public abstract int size();

  /**
   * Test if is symmetric.
   * @return true iff symmetric.
   */
  public abstract boolean isSymmetric();

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < size(); i++) {
      sb.append("[").append(i).append("]");
      final int hi = isSymmetric() ? i : size() - 1;
      for (int j = 0; j <= hi; j++) {
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
