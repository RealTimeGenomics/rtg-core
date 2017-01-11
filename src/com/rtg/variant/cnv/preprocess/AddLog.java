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

import com.rtg.util.MathUtils;

/**
 * Applies log base 2 of a column
 */
public class AddLog implements DatasetProcessor {

  protected final int mCol;
  private final String mColName;

  /**
   * Constructor
   * @param col index of the column to operate on
   * @param colName the name of the new column
   */
  public AddLog(int col, String colName) {
    mCol = col;
    mColName = colName;
  }

  /**
   * Constructor
   * @param col index of the column to operate on
   */
  public AddLog(int col) {
    this(col, null);
  }

  @Override
  public void process(RegionDataset dataset) {
    final NumericColumn in = dataset.asNumeric(mCol);
    final NumericColumn out = dataset.addColumn(new NumericColumn(mColName == null ? "log2(" + dataset.columnName(mCol) + ")" : mColName));
    final double[] values = new double[in.size()];
    for (int i = 0; i < in.size(); ++i) {
      values[i] = Math.log(in.get(i)) / MathUtils.LOG_2;
    }
    out.set(values);
  }
}
