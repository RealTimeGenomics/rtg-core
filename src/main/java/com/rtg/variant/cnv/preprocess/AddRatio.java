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
