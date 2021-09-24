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
