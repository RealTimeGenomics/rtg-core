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

import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 * Utilities for analysing somatic vcf records.
 */
public final class SomaticRecordUtils {
  private SomaticRecordUtils() { }

  /**
   * @param record the record to inspect
   * @return true if the record represents a gain of reference
   */
  public static boolean isGainOfReference(VcfRecord record) {
    final int[] cancerGts = VcfUtils.splitGt(record.getSampleString(AbstractSomaticCaller.CANCER, VcfFormatField.GT.name()));
    return Arrays.stream(cancerGts).allMatch(cancerGt -> cancerGt == 0);
  }

  /**
   * @param record the record to inspect
   * @return true if the record represents a loss of heterozygosity
   */
  public static boolean isLossOfHeterozygosity(VcfRecord record) {
    return lossOfHeterozygosity(record) > 0;
  }

  /**
   * Determine the loss of heterozygosity value
   * &gt; 0 implies loss of heterozygosity
   * 0 no change between normal an cancer
   * &lt; 0 implies
   *
   * @param record vcf record to inspect
   * @return a loss of heterozygosity indicator
   */
  public static double lossOfHeterozygosity(VcfRecord record) {
    final Set<Integer> normalGts = Arrays.stream(VcfUtils.splitGt(record.getSampleString(AbstractSomaticCaller.NORMAL, VcfFormatField.GT.name())))
      .boxed()
      .collect(Collectors.toSet());
    final Set<Integer> cancerGts = Arrays.stream(VcfUtils.splitGt(record.getSampleString(AbstractSomaticCaller.CANCER, VcfFormatField.GT.name())))
      .boxed()
      .collect(Collectors.toSet());
    if (normalGts.equals(cancerGts)) {
      return 0.0;
    }
    return cancerGts.size() < normalGts.size() ? 1.0 : -1.0;
  }

}
