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

package com.rtg.variant.avr;

import java.io.DataOutputStream;
import java.io.IOException;

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.annotation.AbstractDerivedAnnotation;
import com.rtg.vcf.annotation.AnnotationDataType;
import com.rtg.vcf.annotation.DerivedAnnotations;
import com.rtg.vcf.header.VcfHeader;

/**
 * Abstract class to use when implementing a derived annotation.
 */
public class DerivedAnnotation implements Annotation {

  private final AbstractDerivedAnnotation mAnnotation;

  /**
   * @param annotation the derived annotation
   */
  public DerivedAnnotation(AbstractDerivedAnnotation annotation) {
    mAnnotation = annotation;
  }

  /**
   * @param annotationName name of the annotation
   */
  public DerivedAnnotation(String annotationName) {
    this(DerivedAnnotations.valueOf(annotationName).getAnnotation());
  }

  @Override
  public String getName() {
    return "DERIVED-" + mAnnotation.getName();
  }

  @Override
  public AnnotationDataType getType() {
    return mAnnotation.getType();
  }

  @Override
  public String checkHeader(VcfHeader header) {
    return mAnnotation.checkHeader(header);
  }

  @Override
  public Object getValue(VcfRecord record, int sampleNumber) {
    return mAnnotation.getValue(record, sampleNumber);
  }

  @Override
  public void save(DataOutputStream dos) throws IOException {
    dos.writeInt(AnnotationLoader.AnnotationType.DERIVED.ordinal());
    dos.writeUTF(mAnnotation.getName());
  }
}
