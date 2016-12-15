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
 * Computes mean normalization of a column
 */
class MeanNormalize implements DatasetProcessor {

  protected final int mCol;

  protected double mMean;

  MeanNormalize(int col) {
    mCol = col;
  }

  @Override
  public void process(RegionDataset dataset) {
    mMean = Double.NaN;
    if (dataset.size() == 0) {
      return;
    }

    final NumericColumn in = dataset.asNumeric(mCol);
    computeMean(dataset);
    final double[] values = new double[in.size()];
    for (int i = 0; i < in.size(); ++i) {
      values[i] = in.get(i) / mMean;
    }
    final NumericColumn out = dataset.addColumn(new NumericColumn(prefix() + "_" + dataset.columnName(mCol)));
    out.set(values);
  }

  protected String prefix() {
    return "meannorm";
  }

  protected void computeMean(RegionDataset dataset) {
    mMean = dataset.mean(mCol);
  }

  public double mean() {
    return mMean;
  }

}
