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

package com.rtg.calibrate;

import java.util.Locale;

import htsjdk.samtools.SAMRecord;

/**
 * Covariate for read arm.  For paired reads this takes the values 0 and 1, for
 * unpaired data it is always 0.
 */
public final class CovariateArm extends CovariateImpl {

  /**
   * Constructor
   */
  public CovariateArm() {
    super(CovariateEnum.ARM.name().toLowerCase(Locale.ROOT), 2);
  }

  @Override
  public int value(final SAMRecord sam, final CalibratorCigarParser parser) {
    if (sam.getReadPairedFlag()) {
      return sam.getSecondOfPairFlag() ? 1 : 0;
    }
    return 0;
  }

  @Override
  public CovariateEnum getType() {
    return CovariateEnum.ARM;
  }

}
