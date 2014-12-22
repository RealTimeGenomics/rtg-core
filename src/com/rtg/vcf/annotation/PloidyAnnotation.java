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

import java.util.ArrayList;

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.VcfHeader;

/**
 * Ploidy derived attribute
 */
public class PloidyAnnotation extends AbstractDerivedFormatAnnotation {

  /**
   * Constructor
   */
  public PloidyAnnotation() {
    super("PD", "Ploidy of sample. 'h'=>haploid, 'd'=>diploid", AnnotationDataType.STRING);
  }

  @Override
  public Object getValue(VcfRecord record, int sampleNumber) {
    final ArrayList<String> sampleValues = record.getFormatAndSample().get(VcfUtils.FORMAT_GENOTYPE);
    Object value;
    if (sampleValues == null || sampleValues.size() < (sampleNumber + 1)) {
      value = null;
    } else {
      final String sValue = sampleValues.get(sampleNumber);
      if (sValue.contains(VcfRecord.MISSING)) {
        value = null;
      } else {
        value = VcfUtils.isDiploid(record, sampleNumber) ? "d" : "h";
      }
    }
    return value;
  }

  @Override
  public String checkHeader(VcfHeader header) {
    return checkHeader(header, null, new String[]{VcfUtils.FORMAT_GENOTYPE});
  }

}
