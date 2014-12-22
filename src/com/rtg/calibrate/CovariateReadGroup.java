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

import com.rtg.sam.ReadGroupUtils;

import net.sf.samtools.SAMRecord;

/**
 * This covariate grows in size as more read groups are encountered.
 * (Calling parse or value with a new read group adds that read group
 * into the set automatically).
 *
 */
public final class CovariateReadGroup extends CovariateStrings {

  /**
   * Constructor
   */
  public CovariateReadGroup() {
    super(CovariateEnum.READGROUP.name().toLowerCase(Locale.ROOT));
  }

  @Override
  public int value(SAMRecord sam, CalibratorCigarParser parser) {
    return parse(ReadGroupUtils.getReadGroup(sam));
  }

  @Override
  public CovariateEnum getType() {
    return CovariateEnum.READGROUP;
  }

}
