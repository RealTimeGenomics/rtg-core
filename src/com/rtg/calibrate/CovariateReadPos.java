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
 * A Covariate variable for the position along the read.
 *
 * Note:
 * This implementation is incorrect, as the position should be with
 * respect to sequencing order, not mapped reference position order.
 * I.e. it should take into account whether the read was mapped in
 * forward / reverse frame, and whether read is the first or second
 * arm (depending on machine type).
 *
 */
public final class CovariateReadPos extends CovariateImpl {

  /**
   * @param maxReadLen expected maximum length of all reads.
   */
  public CovariateReadPos(int maxReadLen) {
    super(CovariateEnum.READPOSITION.name().toLowerCase(Locale.ROOT), maxReadLen);
  }

  @Override
  public String name() {
    return super.name() + ":" + size();
  }

  @Override
  public int value(SAMRecord sam, CalibratorCigarParser parser) {
    final int readPos = parser.getReadPosition();
    if (readPos >= newSize()) {
      setNewSize(readPos + 1);
    }
    return readPos;
  }

  @Override
  public CovariateEnum getType() {
    return CovariateEnum.READPOSITION;
  }
}
