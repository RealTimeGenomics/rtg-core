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
package com.rtg.variant.bayes.complex;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class SingleCounts extends IntegralAbstract {

  private int mCount = 0;

  private double mCorrection = 0.0;

  /**
   * @param corr correction
   */
  public void increment(final double corr) {
    assert corr >= 0.0 && !Double.isInfinite(corr) && !Double.isNaN(corr);
    ++mCount;
    mCorrection += corr;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(mCount).append(":").append(Utils.realFormat(mCorrection, 3));
  }

  /**
   * @return count
   */
  public int count() {
    return mCount;
  }

  /**
   * @return correction
   */
  public double correction() {
    return mCorrection;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mCount >= 0);
    Exam.assertTrue(mCorrection >= 0.0 && !Double.isInfinite(mCorrection) && !Double.isNaN(mCorrection));
    return true;
  }

}
