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
 * Computes ratio of two columns
 */
public class AddRatio implements DatasetProcessor {

  protected final int mNumCol;
  protected final int mDenomCol;
  private final String mColName;

  /**
   * Constructor
   * @param num index of the column to  use as numerator
   * @param denom index of the column to use as denominator
   * @param colName the name of the new column
   */
  public AddRatio(int num, int denom, String colName) {
    mNumCol = num;
    mDenomCol = denom;
    mColName = colName;
  }

  /**
   * Constructor
   * @param num index of the column to  use as numerator
   * @param denom index of the column to use as denominator
   */
  public AddRatio(int num, int denom) {
    this(num, denom, null);
  }

  @Override
  public void process(RegionDataset dataset) {
    final NumericColumn num = dataset.asNumeric(mNumCol);
    final NumericColumn denom = dataset.asNumeric(mDenomCol);
    assert num.size() == denom.size();
    final NumericColumn out = dataset.addColumn(new NumericColumn(mColName == null ? dataset.columnName(mNumCol) + "/" + dataset.columnName(mDenomCol) : mColName));
    final double[] values = new double[num.size()];
    for (int i = 0; i < num.size(); ++i) {
      values[i] = num.get(i) / denom.get(i);
    }
    out.set(values);
  }
}
