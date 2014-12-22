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
     * Return the value for the given dimension in the SVD.
     * @param dimension to get value for
     * @return value in given dimension
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
   * @param dimensions dimension of SVD matrix
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
   * Set a name for dimension i of the matrix
   * @param i dimension
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
    for (int y = 0; y < mDim; y++) {
      for (int x = 0; x < mDim; x++) {
        if (get(x, y) != get(y, x)) {
          throw new IllegalStateException("Value at (" + x + "," + y + ") = " + get(x, y) + " is not reflected: " + get(y, x));
        }
      }
    }
  }

  private SimilaritySvd removeConstants() {
    final ArrayList<Integer> skiplist = new ArrayList<>();
    for (int j = 0; j < mDim; j++) {
      boolean same = true;
      final long first = get(0, j);
      for (int i = 0; i < mDim; i++) {
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
    for (int j = 0; j < mDim; j++) {
      if (!skiplist.contains(j)) {
        int redi = 0;
        for (int i = 0; i < mDim; i++) {
          if (!skiplist.contains(i)) {
            reduced.put(redi, redj, get(i, j));
            redi++;
          }
        }
        reduced.putName(redj, getName(j));
        redj++;
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
      // to few dimensions to process
      return;
    }

    double t = 0.0;
    final double[] values = new double[n * n];
    for (int j = 0; j < n; j++) {
      for (int i = 0; i <= j; i++) {
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
    for (int j = 0; j < n; j++) {
      for (int i = 0; i <= j; i++) {
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
    for (int j = 0; j < n; j++) {
      final SvdRow sr = new SvdRow(dimensions, mat.getName(j));
      for (int i = 0; i < dimensions; i++) {
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
   * Return the value for row, dimension from SVD decomposition.
   * @param row row in decomposition
   * @param dimension dimension in decomposition
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
    for (int j = 0; j < mDim; j++) {
      for (int i = 0; i < mDim; i++) {
        sb.append(get(i, j)).append(StringUtils.TAB);
      }
      sb.append(getName(j)).append(StringUtils.LS);
    }
    return sb.toString();
  }
}
