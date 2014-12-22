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

package com.rtg.vcf.annotation;

import java.util.List;

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.VcfHeader;

/**
 * Derived attribute giving GQ divided by DP for a sample.
 */
public class GenotypeQualityOverDepthAnnotation extends AbstractDerivedFormatAnnotation {

  /**
   * Constructor.
   */
  public GenotypeQualityOverDepthAnnotation() {
    super("GQD", "GQ / DP for a single sample", AnnotationDataType.DOUBLE);
  }

  @Override
  public Object getValue(VcfRecord record, int sampleNumber) {
    final List<String> dps = record.getFormatAndSample().get(VcfUtils.FORMAT_SAMPLE_DEPTH);
    final List<String> gqs = record.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE_QUALITY);
    if (dps != null && gqs != null && dps.size() > sampleNumber && dps.size() == gqs.size()) {
      final String dpVal = dps.get(sampleNumber);
      final String gqVal = gqs.get(sampleNumber);
      if (!VcfRecord.MISSING.equals(dpVal) && !VcfRecord.MISSING.equals(gqVal)) {
        final int dp = Integer.parseInt(dpVal);
        if (dp <= 0) {
          return Double.POSITIVE_INFINITY;
        }
        return Double.valueOf(gqVal) / dp;
      }
    }
    return null;
  }

  @Override
  public String checkHeader(VcfHeader header) {
    return checkHeader(header, null, new String[]{VcfUtils.FORMAT_SAMPLE_DEPTH, VcfUtils.FORMAT_GENOTYPE_QUALITY});
  }
}
