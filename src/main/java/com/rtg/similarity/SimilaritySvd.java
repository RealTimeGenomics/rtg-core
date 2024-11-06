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

package com.rtg.similarity;

import java.util.ArrayList;

import com.rtg.util.StringUtils;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

/**
 * Computes single value decomposition for a matrix.
 *
 */
public class SimilaritySvd {
  /**
   * A row in the principle component analysis output.
   */
  public static final class SvdRow {
    private final double[] mValues;
    private final String mName;

    private SvdRow(int dimensions, String name) {
      assert dimensions >= 1;
      assert name != null;
      mValues = new double[dimensions];
      mName = name;
    }

    public int getDimensions() {
      return mValues.length;
    }

    /**
     * Return the value for the given length in the SVD.
     * @param dimension to get value for
     * @return value in given length
     */
    public double getValue(int dimension) {
      return mValues[dimension];
    }

    private void setValue(int dimension, double value) {
      mValues[dimension] = value;
    }

    /**
     * Get name associated with row in SVD.
     * @return name of row
     */
    public String getName() {
      return mName;
    }
  }

  private final int mDim;
  private final long[] mValues;
  private final String[] mNames;

  private SvdRow[] mRows = null;

  /**
   * Create a similarity SVD with the given number of dimensions
   * @param dimensions length of SVD matrix
   */
  public SimilaritySvd(int dimensions) {
    mDim = dimensions;
    mValues = new long[dimensions * dimensions];
    mNames = new String[dimensions];
  }

  /**
   * Put a value in position i, j of the matrix
   * @param i position in x
   * @param j position in y
   * @param value value to set
   */
  public void put(int i, int j, long value) {
    mValues[j * mDim + i] = value;
  }

  private long get(int i, int j) {
    return mValues[j * mDim + i];
  }

  /**
   * Set a name for length i of the matrix
   * @param i length
   * @param name name to set
   */
  public void putName(int i, String name) {
    mNames[i] = name;
  }

  private String getName(int i) {
    return mNames[i];
  }

  private int getDimension() {
    return mDim;
  }

  private void checkConsistent() {
    for (int y = 0; y < mDim; ++y) {
      for (int x = 0; x < mDim; ++x) {
        if (get(x, y) != get(y, x)) {
          throw new IllegalStateException("Value at (" + x + "," + y + ") = " + get(x, y) + " is not reflected: " + get(y, x));
        }
      }
    }
  }

  private SimilaritySvd removeConstants() {
    final ArrayList<Integer> skiplist = new ArrayList<>();
    for (int j = 0; j < mDim; ++j) {
      boolean same = true;
      final long first = get(0, j);
      for (int i = 0; i < mDim; ++i) {
        if (get(i, j) != first) {
          same = false;
        }
      }
      if (same) {
        //System.err.println("Row " + (j + 1) + " '" + getName(j) + "' is constant. Values are all " + first + ". Removing row/column.");
        skiplist.add(j);
      }
    }

    if (skiplist.size() == 0) {
      return this;
    }

    final SimilaritySvd reduced = new SimilaritySvd(mDim - skiplist.size());

    int redj = 0;
    for (int j = 0; j < mDim; ++j) {
      if (!skiplist.contains(j)) {
        int redi = 0;
        for (int i = 0; i < mDim; ++i) {
          if (!skiplist.contains(i)) {
            reduced.put(redi, redj, get(i, j));
            ++redi;
          }
        }
        reduced.putName(redj, getName(j));
        ++redj;
      }
    }
    return reduced;
  }

  /**
   * Decomposes the given matrix down to the specified dimensions.
   *
   * @param dimensions number of dimensions to reduce to
   */
  public void decompose(int dimensions) {
    checkConsistent();

    final SimilaritySvd mat = removeConstants();

    // find average of normalized values
    final int n = mat.getDimension();

    if (n <= dimensions) {
      // too few dimensions to process
      return;
    }

    double t = 0.0;
    final double[] values = new double[n * n];
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i <= j; ++i) {
        final double d = mat.get(i, j) / (double) (mat.get(j, j) + mat.get(i, i) - mat.get(i, j));
        values[j * n + i] = d;
        if (d < 1.0) {
          t += d;
        }
      }
    }
    t *= 2.0;
    final double ave = t / (n * n - n);

    //System.err.println(mat);

    // create N-by-N matrix that doesn't have full rank
    // normalized values - average values
    final Matrix matA = new Matrix(n, n);
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i <= j; ++i) {
        final double d = values[j * n + i] - ave;
        matA.set(i, j, d);
        matA.set(j, i, d);
      }
    }

    // compute the singular value decomposition
    //System.err.println("calculating SVD " + getDimension() + " " + mat.getDimension());
    final SingularValueDecomposition s = matA.svd();
    final Matrix matUk = s.getU().getMatrix(0, n - 1, 0, dimensions);

    final SvdRow[] res = new SvdRow[n];
    for (int j = 0; j < n; ++j) {
      final SvdRow sr = new SvdRow(dimensions, mat.getName(j));
      for (int i = 0; i < dimensions; ++i) {
        sr.setValue(i, matUk.get(j, i));
      }
      res[j] = sr;
    }
    mRows = res;
  }

  /**
   * Return number of rows in SVD decomposition
   * @return number of rows, or -1 if not decomposed.
   */
  public int getSvdRowsLength() {
    return mRows == null ? -1 : mRows.length;
  }

  /**
   * Return dimensions in SVD decomposition
   * @return number of dimensions, or -1 if not decomposed.
   */
  public int getSvdDimension() {
    return getSvdRowsLength() <= 0 ? -1 : mRows[0].getDimensions();
  }

  /**
   * Return the value for row, length from SVD decomposition.
   * @param row row in decomposition
   * @param dimension length in decomposition
   * @return PCA value
   */
  public double getSvdValue(int row, int dimension) {
    if (getSvdRowsLength() < 0) {
      throw new IllegalStateException("SVD decomposition not computed.");
    }
    return mRows[row].getValue(dimension);
  }

  /**
   * Returns the name of the row in the SVD decomposition.
   * @param row row in decomposition
   * @return row name
   */
  public String getSvdName(int row) {
    if (getSvdRowsLength() < 0) {
      throw new IllegalStateException("SVD decomposition not computed.");
    }
    return mRows[row].getName();
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (int j = 0; j < mDim; ++j) {
      for (int i = 0; i < mDim; ++i) {
        sb.append(get(i, j)).append(StringUtils.TAB);
      }
      sb.append(getName(j)).append(StringUtils.LS);
    }
    return sb.toString();
  }
}
