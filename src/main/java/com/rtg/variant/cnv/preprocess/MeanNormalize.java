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
    final NumericColumn out = dataset.addColumn(new NumericColumn(prefix() + "(" + dataset.columnName(mCol) + ")"));
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
