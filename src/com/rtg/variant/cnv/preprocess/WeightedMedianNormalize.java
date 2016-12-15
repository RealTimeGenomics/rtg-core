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
 * Computes median normalization of a column, incorporating region lengths into median calculation
 */
public class WeightedMedianNormalize extends MedianNormalize {

  /**
   * Constructor
   * @param col index of the column to operate on
   */
  public WeightedMedianNormalize(int col) {
    super(col);
  }

  @Override
  protected void computeMedian(RegionDataset dataset) {
    mMedian = dataset.weightedMedian(mCol);
  }

  @Override
  protected String prefix() {
    return "weightedmediannorm";
  }

}
