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

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Applies GC median normalization to a column.
 * The dataset must already have a column containing the GC percent of the region.
 */
public class GcNormalize implements DatasetProcessor {

  private final int mBins;
  private final double mInterval;
  private final int mCol;
  private final String mColName;

  /**
   * Constructor
   * @param col index of the column to operate on
   * @param bins number of bins to use
   * @param colName name to use for the new column
   */
  public GcNormalize(int col, int bins, String colName) {
    mCol = col;
    mBins = bins;
    mInterval = 1.0 / mBins;
    mColName = colName;
  }

  /**
   * Constructor
   * @param col index of the column to operate on
   * @param bins number of bins to use
   */
  public GcNormalize(int col, int bins) {
    this(col, bins, null);
  }

  private int bin(double gc) {
    final int bin = (int) (gc / mInterval);
    return bin == mBins ? mBins - 1 : bin;
  }

  private double binMid(int bin) {
    return bin * mInterval + mInterval / 2;
  }

  @Override
  public void process(RegionDataset dataset) {
    final int gcCol = dataset.columnId(AddGc.PCT_GC_NAME);
    if (gcCol < 0) {
      throw new NoTalkbackSlimException("Dataset does not contain " + AddGc.PCT_GC_NAME + " column. " + dataset.getColumnNames());
    }
    final NumericColumn gc = dataset.asNumeric(gcCol);
    final NumericColumn in = dataset.asNumeric(mCol);
    final NumericColumn out = dataset.addColumn(new NumericColumn(mColName == null ? "gcnorm(" + dataset.columnName(mCol) + ")" : mColName));

    final double median = in.median();

    final NumericColumn[] binned = new NumericColumn[mBins];
    for (int i = 0; i < binned.length; ++i) {
      binned[i] = new NumericColumn(in.getName());
    }
    for (int i = 0; i < in.size(); ++i) {
      binned[bin(gc.get(i))].add(in.get(i));
    }
    final double[] norms = new double[binned.length];
    for (int i = 0; i < binned.length; ++i) {
      final NumericColumn c = binned[i];
      if (c.size() == 0) {
        continue;
      }
      final double binmedian = c.median();
      norms[i] = binmedian > 0 ? median / binmedian : 1.0;
      Diagnostic.developerLog(String.format("GC Bin %d n=%d %%gc=%.2f median=%.4f norm=%.4f", i, c.size(), binMid(i), binmedian, norms[i]));
    }

    final double[] values = new double[dataset.size()];
    for (int i = 0; i < values.length; ++i) {
      values[i] = in.get(i) * norms[bin(gc.get(i))];
    }
    out.set(values);
  }
}
