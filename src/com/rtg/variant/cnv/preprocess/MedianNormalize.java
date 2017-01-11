/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.cnv.preprocess;

/**
 * Computes median normalization of a column
 */
public class MedianNormalize implements DatasetProcessor {

  /** Column containing input data */
  protected final int mCol;
  private final String mColName;

  protected double mMedian;

  /**
   * Constructor
   * @param col index of the column to operate on
   * @param colName the name of the new column
   */
  public MedianNormalize(int col, String colName) {
    mCol = col;
    mColName = colName;
  }

  /**
   * Constructor
   * @param col index of the column to operate on
   */
  public MedianNormalize(int col) {
    this(col, null);
  }

  @Override
  public void process(RegionDataset dataset) {
    mMedian = Double.NaN;
    if (dataset.size() == 0) {
      return;
    }
    final NumericColumn in = dataset.asNumeric(mCol);
    computeMedian(dataset);
    final double[] values = new double[in.size()];
    for (int i = 0; i < in.size(); ++i) {
      values[i] = in.get(i) / mMedian;
    }
    final NumericColumn out = dataset.addColumn(new NumericColumn(mColName == null ? prefix() + "(" + dataset.columnName(mCol) + ")" : mColName));
    out.set(values);
  }

  protected void computeMedian(RegionDataset dataset) {
    mMedian = dataset.median(mCol);
  }

  protected String prefix() {
    return "mediannorm";
  }

  protected double median() {
    return mMedian;
  }
}
