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
package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.VcfFilter;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

/**
 * A filter serving two functions in somatic calling. It ensures that all records are passed
 * into the SomaticStatistics for accumulation of contamination estimates and controls
 * whether or not non-somatic variants are output.
 */
public class SomaticFilter implements VcfFilter {

  private final SomaticStatistics mStatistics;
  private final boolean mSomaticOnly;
  private VcfHeader mVcfHeader;
  private final boolean mLossOfHeterozygosity;
  private final boolean mGainOfReference;

  /**
   * Construct a new filter which counts variants for the purposes of estimating contamination in
   * the somatic caller and optionally filter records which are not somatic.
   * @param statistics statistics object for computing contamination estimates
   * @param somaticOnly true iff only somatic variants are to be retained
   * @param lossOfHeterozygosity true iff loss of heterozygosity calls should be retained
   * @param gainOfReference true iff loss of heterozygosity calls should be retained
   */
  public SomaticFilter(final SomaticStatistics statistics, final boolean somaticOnly, final boolean lossOfHeterozygosity, final boolean gainOfReference) {
    mStatistics = statistics;
    mSomaticOnly = somaticOnly;
    mLossOfHeterozygosity = lossOfHeterozygosity;
    mGainOfReference = gainOfReference;
  }

  @Override
  public void setHeader(VcfHeader vcfHeader) {
    mVcfHeader = vcfHeader;
  }

  @Override
  public boolean accept(VcfRecord record) {
    mStatistics.countVariant(mVcfHeader, record);
    final Integer somaticStatus = record.getSampleInteger(AbstractSomaticCaller.CANCER, VcfFormatField.SS.name());
    final boolean somatic = somaticStatus != null && somaticStatus == 2;
    if (mSomaticOnly && !somatic) {
      return false;

    }

    if (somatic && !mLossOfHeterozygosity && SomaticRecordUtils.isLossOfHeterozygosity(record)) {
      return false;
    }
    if (somatic && !mGainOfReference && SomaticRecordUtils.isGainOfReference(record)) {
      return false;
    }
    return true;
  }

}
