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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.bed.BedRangeLoader;
import com.rtg.bed.BedRecord;

/**
 * Read a BED file containing site specific prior expressed as a probability in column 4.
 */
public class PriorBedRangeLoader extends BedRangeLoader<Double> {

  /**
   * Construct a new <code>BedRangeLoader</code> for a BED giving priors.
   */
  public PriorBedRangeLoader() {
    super(1);
  }

  @Override
  protected Double getMeta(BedRecord rec) {
    if (rec.getAnnotations() != null && rec.getAnnotations().length > 0) {
      final double p = Double.parseDouble(rec.getAnnotations()[0]);
      assert p >= 0 && p <= 1;
      return p;
    }
    return Double.NaN;
  }
}
