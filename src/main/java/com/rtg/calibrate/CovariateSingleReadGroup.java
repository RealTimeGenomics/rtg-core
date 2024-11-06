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

package com.rtg.calibrate;

import java.util.Locale;

import htsjdk.samtools.SAMRecord;

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
