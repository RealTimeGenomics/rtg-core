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

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.MathUtils;

/**
 * Applies log base 2 of a column
 */
public class AddLog implements DatasetProcessor {

  private static final double MIN_LOG_RATIO = GlobalFlags.getDoubleValue(CoreGlobalFlags.SEGMENT_MIN_LOG_RATIO);

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
      values[i] = Math.max(Math.log(in.get(i)) / MathUtils.LOG_2, MIN_LOG_RATIO);
    }
    out.set(values);
  }
}
