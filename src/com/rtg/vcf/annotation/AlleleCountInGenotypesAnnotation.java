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
 * Implementation of the VCF format specified AC field.
 * Allele count in genotypes, for each alternative allele, in the same order as listed.
 */
public class AlleleCountInGenotypesAnnotation extends AbstractDerivedAnnotation {

  /**
   * Constructor
   */
  public AlleleCountInGenotypesAnnotation() {
    super("AC", "Allele count in genotypes, for each alternative allele, in the same order as listed", AnnotationDataType.INTEGER);
  }

  @Override
  public Object getValue(VcfRecord record, int sampleNumber) {
    final List<String> gtList = record.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE);
    if (gtList == null || record.getAltCalls().size() == 0) {
      return null;
    }
    final int[] ret = new int[record.getAltCalls().size()];
    for (final String gtStr : gtList) {
      final int[] gts = VcfUtils.splitGt(gtStr);
      for (final int gt : gts) {
        if (gt > 0) {
          ret[gt - 1]++;
        }
      }
    }
    return ret;
  }

  @Override
  public String checkHeader(VcfHeader header) {
    return checkHeader(header, null, new String[]{VcfUtils.FORMAT_GENOTYPE});
  }

}
