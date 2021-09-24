/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.sv;

import com.rtg.util.Utils;


/**
 * Constant.
 */
public class DistributionArray extends Distribution {

  private final double[] mValues;

  /**
   * @param lo low index of the distribution (inclusive).
   * @param values used in underlying array.
   */
  public DistributionArray(final int lo, final double[] values) {
    super(lo, lo + values.length);
    mValues = values.clone();
    globalIntegrity();
  }

  @Override
  protected double getValue(final int index) {
    return mValues[index - lo()];
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(Utils.realFormat(mValues, 3));
  }

  @Override
  public boolean integrity() {
    super.integrity();
    return true;
  }

}
