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

package com.rtg.vcf;

import static com.rtg.vcf.VcfFilterStatistics.Stat;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.vcf.header.VcfHeader;

/**
 * Abstract out process of checking a condition and recording stats about it
 */
@TestClass("com.rtg.vcf.VcfFilterCliTest")
public abstract class AbstractVcfFilter implements VcfFilter {
  final VcfFilterStatistics mStatistics;
  final Stat mStat;

  AbstractVcfFilter(VcfFilterStatistics stats, Stat stat) {
    mStatistics = stats;
    mStat = stat;
  }

  @Override
  public boolean accept(VcfRecord record) {
    if (!acceptCondition(record)) {
      mStatistics.increment(mStat);
      return false;
    }
    return true;
  }

  @Override
  public void setHeader(VcfHeader header) {
  }

  /**
   * Condition to filter on
   * @param record the VCF record to filter
   * @return false if the record is unacceptable and should be filtered
   */
  abstract boolean acceptCondition(VcfRecord record);

  /**
   * Filter on QUAL field
   */
  public static class QualFilter extends AbstractVcfFilter {
    final double mMinQuality;
    final double mMaxQuality;
    QualFilter(VcfFilterStatistics stats, double minQuality, double maxQuality) {
      super(stats, Stat.QUALITY_FILTERED_COUNT);
      mMinQuality = minQuality;
      mMaxQuality = maxQuality;
    }
    @Override
    boolean acceptCondition(VcfRecord record) {
      // QUAL filtering
      final String qualityStr = record.getQuality();
      if (!VcfRecord.MISSING.equals(qualityStr)) {
        final double quality = Double.parseDouble(qualityStr);
        if (quality < mMinQuality || quality > mMaxQuality) {
          return false;
        }
      }
      return true;
    }
  }

}
