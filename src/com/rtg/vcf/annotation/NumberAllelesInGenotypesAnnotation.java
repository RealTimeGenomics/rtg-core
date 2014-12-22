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
 * Implementation of the AN info field as defined in the VCF specification.
 * Total number of alleles in called genotypes.
 */
public class NumberAllelesInGenotypesAnnotation extends AbstractDerivedAnnotation {

  /**
   * Constructor
   */
  public NumberAllelesInGenotypesAnnotation() {
    super("AN", "Total number of alleles in called genotypes", AnnotationDataType.INTEGER);
  }

  @Override
  public Object getValue(VcfRecord record, int sampleNumber) {
    final List<String> gtList = record.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE);
    if (gtList == null) {
      return null;
    }
    int count = 0;
    for (final String gtStr : gtList) {
      final int[] gts = VcfUtils.splitGt(gtStr);
      for (final int gt : gts) {
        if (gt >= 0) {
          count++;
        }
      }
    }
    return count;
  }

  @Override
  public String checkHeader(VcfHeader header) {
    return checkHeader(header, null, new String[]{VcfUtils.FORMAT_GENOTYPE});
  }

}
