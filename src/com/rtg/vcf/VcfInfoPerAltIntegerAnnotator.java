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

import com.rtg.vcf.annotation.AbstractDerivedAnnotation;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 */
public class VcfInfoPerAltIntegerAnnotator implements VcfAnnotator {

  final AbstractDerivedAnnotation mAnnotation;

  /**
   * Create an INFO annotation that outputs an integer value for each alt allele in the variant.
   * @param annotation the annotation to use.
   */
  public VcfInfoPerAltIntegerAnnotator(AbstractDerivedAnnotation annotation) {
    assert annotation != null && annotation.getType().getClassType() == Integer.class;
    mAnnotation = annotation;
  }

  @Override
  public void updateHeader(VcfHeader header) {
    header.ensureContains(new InfoField(mAnnotation.getName(), MetaType.INTEGER, VcfNumber.ALTS, mAnnotation.getDescription()));
  }

  @Override
  public void annotate(VcfRecord rec) {
    final int[] vals = (int[]) mAnnotation.getValue(rec, -1);
    if (vals != null) {
      final String[] vcfVals = new String[vals.length];
      for (int i = 0 ; i < vals.length; i++) {
        vcfVals[i] = Integer.toString(vals[i]);
      }
      rec.setInfo(mAnnotation.getName(), vcfVals);
    }
  }

}
