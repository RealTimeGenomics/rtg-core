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
import com.rtg.vcf.header.VcfHeader;

/**
 * Annotation that counts the number of alternative alleles in a VCF record.
 */
public class NumberOfAltAllelesAnnotation extends AbstractDerivedAnnotation {

  /**
   * Constructor.
   */
  public NumberOfAltAllelesAnnotation() {
    super("NAA", "Number of alternative alleles", AnnotationDataType.INTEGER);
  }

  @Override
  public Integer getValue(VcfRecord record, int sampleNumber) {
    final List<String> alleles = record.getAltCalls();
    int missingCount = 0;
    for (String allele : alleles) {
      if (VcfRecord.MISSING.equals(allele)) {
        missingCount++;
      }
    }
    return record.getAltCalls().size() - missingCount;
  }

  @Override
  public String checkHeader(VcfHeader header) {
    // Alleles column is not optional
    return null;
  }

}
