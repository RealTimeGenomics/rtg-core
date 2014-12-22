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

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.VcfHeader;

/**
 * Annotation for the length of the longest allele across all samples in VCF record.
 */
public class LongestAlleleAnnotation extends AbstractDerivedAnnotation {

  /**
   * Constructor
   */
  public LongestAlleleAnnotation() {
    super("LAL", "Length of longest allele", AnnotationDataType.INTEGER);
  }

  @Override
  public Object getValue(VcfRecord record, int sampleNumber) {
    int length = record.getRefCall().length();
    for (final String allele : record.getAltCalls()) {
      length = Math.max(length, allele.length());
    }
    return length;
  }

  @Override
  public String checkHeader(VcfHeader header) {
    return null; //Alleles column is not optional
  }

}
