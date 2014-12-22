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
 * Derived attribute giving QUAL divided by DP.
 */
public class QualOverDepthAnnotation extends AbstractDerivedAnnotation {

  /**
   * Constructor.
   */
  public QualOverDepthAnnotation() {
    super("QD", "QUAL / DP", AnnotationDataType.DOUBLE);
  }

  @Override
  public Double getValue(VcfRecord record, int sampleNumber) {
    final String squal = record.getQuality();
    if (!VcfRecord.MISSING.equals(squal)) {
      if (record.getInfo().containsKey(VcfUtils.INFO_COMBINED_DEPTH)) {
        // get dp from info DP field
        final String sdp = record.getInfo().get(VcfUtils.INFO_COMBINED_DEPTH).get(0);
        if (sdp != null && !VcfRecord.MISSING.equals(sdp)) {
          final int dp = Integer.parseInt(sdp);
          if (dp <= 0) {
            return Double.POSITIVE_INFINITY;
          }
          return Double.parseDouble(squal) / dp;
        }
      } else {
        final List<String> dps = record.getFormatAndSample().get(VcfUtils.FORMAT_SAMPLE_DEPTH);
        if (dps != null && dps.size() != 0) {
          // get dp from sum of DP in samples, as DP may not exist in INFO for single sample VCF
          int dp = 0;
          for (String sdp : dps) {
            if (!VcfRecord.MISSING.equals(sdp)) {
              dp += Integer.parseInt(sdp);
            }
          }
          if (dp <= 0) {
            return Double.POSITIVE_INFINITY;
          }
          return Double.parseDouble(squal) / dp;
        }
      }
    }
    return null;
  }

  @Override
  public String checkHeader(VcfHeader header) {
    // QUAL column is not optional
    final String noCombinedDP = checkHeader(header, new String[]{VcfUtils.INFO_COMBINED_DEPTH}, null);
    final String noSampleDP = checkHeader(header, null, new String[]{VcfUtils.FORMAT_SAMPLE_DEPTH});
    // Only need one or the other, not both, but if neither is present, ask for combined depth column
    // as that is what this is technically supposed to run from
    return noCombinedDP != null && noSampleDP != null ? noCombinedDP : null;
  }
}
