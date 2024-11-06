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

import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 * Utilities for analysing somatic VCF records.
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
   * @param record VCF record to inspect
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
