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

import com.rtg.sam.BadSuperCigarException;

import htsjdk.samtools.SAMRecord;

/**
 * A Covariate variable for the quality of the read base.
 *
 */
public final class CovariateBaseQuality extends CovariateImpl {

  private static final int MAX_QUAL = 64;

  /**
   * Create a Covariate variable for base qualities.
   */
  public CovariateBaseQuality() {
    super(CovariateEnum.BASEQUALITY.name().toLowerCase(Locale.ROOT), MAX_QUAL);
  }

  @Override
  public int value(SAMRecord sam, CalibratorCigarParser parser) throws BadSuperCigarException {
    final int ret = parser.getCurrentQuality();
    if (ret >= newSize()) {
      setNewSize(ret + 1);
    }
    return ret;
  }

  @Override
  public CovariateEnum getType() {
    return CovariateEnum.BASEQUALITY;
  }
}
