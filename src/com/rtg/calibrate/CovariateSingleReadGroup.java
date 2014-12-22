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

import net.sf.samtools.SAMRecord;

/**
 */
public final class CovariateSingleReadGroup extends CovariateImpl {

  private final String mReadGroup;

  /**
   * Covariate for read group in the case where we know in advance there is only one read group.
   * @param readGroup the name of the read group.
   */
  public CovariateSingleReadGroup(final String readGroup) {
    super(CovariateEnum.READGROUP.name().toLowerCase(Locale.ROOT), 1);
    mReadGroup = readGroup;
  }

  @Override
  public String valueString(final int v) {
    return mReadGroup;
  }

  @Override
  public int value(SAMRecord sam, CalibratorCigarParser parser) {
    return 0;
  }

  @Override
  public int parse(final String v) {
    if (mReadGroup.equals(v)) {
      return 0;
    }
    throw new UnsupportedOperationException();
  }

  @Override
  public CovariateEnum getType() {
    return CovariateEnum.READGROUP;
  }

}
