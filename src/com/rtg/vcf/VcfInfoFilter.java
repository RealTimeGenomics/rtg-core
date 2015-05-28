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

import java.util.ArrayList;

import com.rtg.vcf.VcfFilterStatistics.Stat;

/**
 */
public abstract class VcfInfoFilter extends AbstractVcfFilter {

  VcfInfoFilter(VcfFilterStatistics stats, Stat stat) {
    super(stats, stat);
  }

  /**
   * Filter on the range of an integer field
   */
  public static class MinMaxIntFilter extends VcfInfoFilter {
    private final String mField;
    private final int mMin;
    private final int mMax;
    MinMaxIntFilter(VcfFilterStatistics stats, Stat stat, int min, int max, String field) {
      super(stats, stat);
      mMin = min;
      mMax = max;
      mField = field;
    }
    @Override
    boolean acceptCondition(VcfRecord record) {
      final ArrayList<String> values = record.getInfo().get(mField);
      if (values != null && values.size() == 1) {
        final Integer value = Integer.valueOf(values.get(0));
        return !(value < mMin || value > mMax);
      }
      return true;
    }
  }
}
