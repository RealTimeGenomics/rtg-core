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
import com.rtg.util.MathUtils;
import com.rtg.variant.PosteriorUtils;
import com.rtg.variant.format.VcfFormatField;

/**
 */
@TestClass("com.rtg.vcf.VcfFilterCliTest")
public abstract class VcfSampleFilter extends AbstractVcfFilter {
  int[] mSamples = null;
  boolean[] mSampleFailed = null;
  VcfSampleFilter(VcfFilterStatistics stats, Stat stat) {
    super(stats, stat);
  }

  void setSamples(int[] samples, boolean[] failedSamples) {
    mSamples = samples;
    mSampleFailed = failedSamples;
  }

  @Override
  boolean acceptCondition(VcfRecord record) {
    boolean result = true;
    if (mSamples != null) {
      for (int i : mSamples) {
        if (!acceptSample(record, i)) {
          result = false;
          if (mSampleFailed != null) {
            mSampleFailed[i] = true;
          } else {
            break;
          }
        }
      }
    }
    return result;
  }

  abstract boolean acceptSample(VcfRecord record, int index);


  /**
   * Filter on the range of a double field
   */
  public static class MinMaxDoubleFilter extends VcfSampleFilter {
    final double mMin;
    final double mMax;
    final String mField;
    MinMaxDoubleFilter(VcfFilterStatistics stats, Stat stat, double min, double max, String field) {
      super(stats, stat);
      mField = field;
      mMin = min;
      mMax = max;
    }
    @Override
    boolean acceptSample(VcfRecord record, int sampleIndex) {
      // check ambiguity ratio
      final Double val = record.getSampleDouble(sampleIndex, mField);
      return val == null || !(val < mMin || val > mMax);
    }
  }

  /**
   * Filter on the range of an integer field
   */
  public static class MinMaxIntFilter extends VcfSampleFilter {
    final int mMin;
    final int mMax;
    final String mField;
    MinMaxIntFilter(VcfFilterStatistics stats, Stat stat, int min, int max, String field) {
      super(stats, stat);
      mField = field;
      mMin = min;
      mMax = max;
    }
    @Override
    boolean acceptSample(VcfRecord record, int sampleIndex) {
      // check ambiguity ratio
      final Integer val = record.getSampleInteger(sampleIndex, mField);
      return val == null || !(val < mMin || val > mMax);
    }
  }

  /**
   * Filter on Genotype Quality
   */
  public static class GqFilter extends VcfSampleFilter {
    final double mMinGq;
    final double mMaxGq;
    final boolean mPosteriorFiltering;
    GqFilter(VcfFilterStatistics stats, double minQuality, double maxQuality, boolean posteriorFiltering) {
      super(stats, Stat.GENOTYPE_QUALITY_POSTERIOR_FILTERED_COUNT);
      mMinGq = minQuality;
      mMaxGq = maxQuality;
      mPosteriorFiltering = posteriorFiltering;
    }
    @Override
    boolean acceptSample(VcfRecord record, int index) {
      // GQ filtering
      final Double gq = record.getSampleDouble(index, VcfUtils.FORMAT_GENOTYPE_QUALITY);
      if (gq != null) {
        if (mPosteriorFiltering) {
          // check posterior
          final double pgq = PosteriorUtils.unphredIfy(gq) / MathUtils.LOG_10;
          if (pgq < mMinGq || pgq > mMaxGq) {
            return false;
          }
        } else {
          // check genotype quality
          if (gq < mMinGq || gq > mMaxGq) {
            return false;
          }
        }
      }
      return true;
    }
  }

  /**
   * Only allow filtering on a single sample.
   * Records will be accepted if that sample has a de Novo call within the score range
   * and it is the only such sample
   */
  public static class DenovoFilter extends VcfSampleFilter {
    final double mMinDenovoScore;
    final double mMaxDenovoScore;
    DenovoFilter(VcfFilterStatistics stats, double minQuality, double maxQuality) {
      super(stats, Stat.DENOVO_SCORE);
      mMinDenovoScore = minQuality;
      mMaxDenovoScore = maxQuality;
    }

    @Override
    boolean acceptCondition(VcfRecord record) {
      assert mSamples.length == 1;
      for (int sampleIndex : mSamples) {
        if (!"Y".equals(record.getSampleString(sampleIndex, VcfFormatField.DN.name()))) {
          if (mSampleFailed != null) {
            mSampleFailed[sampleIndex] = true;
          }
          return false;
        }
      }
      boolean result = false;
      for (int sampleIndex : mSamples) {
        final Double dnp = record.getSampleDouble(sampleIndex, VcfFormatField.DNP.name());
        if (dnp != null) {
          if (dnp >= mMinDenovoScore && dnp <= mMaxDenovoScore) {
            result = true;
            for (int i = 0; i < record.getNumberOfSamples(); i++) {
              final Double otherDnp = record.getSampleDouble(i, VcfFormatField.DNP.name());
              final String otherDn = record.getSampleString(i, VcfFormatField.DN.name());
              if (i != sampleIndex && "Y".equals(otherDn) && otherDnp != null && otherDnp >= mMinDenovoScore && otherDnp <= mMaxDenovoScore) {
                result = false;
              }
            }
            return result;
          } else {
            if (mSampleFailed != null) {
              mSampleFailed[sampleIndex] = true;
            }
            return false;
          }

        }
      }
      return result;
    }
    @Override
    boolean acceptSample(VcfRecord record, int index) {
      throw new IllegalArgumentException("De novo filter is a bit weird don't call this");
    }

  }
}
