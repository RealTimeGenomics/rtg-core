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

import java.io.IOException;

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

  /**
   * Constructor
   * @param col index of the column to operate on
   * @param bins number of bins to use
   */
  public GcNormalize(int col, int bins) {
    mCol = col;
    mBins = bins;
    mInterval = 1.0 / mBins;
  }

  private int bin(double gc) {
    final int bin = (int) (gc / mInterval);
    return bin == mBins ? mBins - 1 : bin;
  }

  private double binMid(int bin) {
    return bin * mInterval + mInterval / 2;
  }

  @Override
  public void process(RegionDataset dataset) throws IOException {
    final int gcCol = dataset.columnId(AddGc.PCT_GC_NAME);
    if (gcCol < 0) {
      throw new NoTalkbackSlimException("Dataset does not contain " + AddGc.PCT_GC_NAME + " column. " + dataset.getColumnNames());
    }
    final NumericColumn gc = dataset.asNumeric(gcCol);
    final NumericColumn in = dataset.asNumeric(mCol);
    final NumericColumn out = dataset.addColumn(new NumericColumn("gcnorm(" + dataset.columnName(mCol) + ")"));

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
