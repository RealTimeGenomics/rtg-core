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

package com.rtg.variant.bayes;

import com.rtg.util.MathUtils;

/**
 */
public final class EntropyBinomial {

  private EntropyBinomial() { }

  /**
   * @param total sum of counts.
   * @param counts of individual categories.
   * @param prior the log of the prior.
   * @return return ln of posterior for random distribution.
   */
  static double noisePosteriorLocal(final int total, final int[] counts, double prior) {
    double sum = -MathUtils.logFactorial(total);
    for (final int c : counts) {
      sum += MathUtils.logFactorial(c);
    }
    assert !Double.isInfinite(sum) && !Double.isNaN(sum);
    return sum + prior;
  }

}
