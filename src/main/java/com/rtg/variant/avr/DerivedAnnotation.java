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
import com.rtg.vcf.annotation.DerivedAnnotations;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;

/**
 * Derived annotations compute a value on the fly based on other existing fields in the record.
 */
public class DerivedAnnotation implements Annotation {

  private final AbstractDerivedAnnotation<?> mAnnotation;
  private final AnnotationDataType mType;

  /**
   * @param annotation the derived annotation
   */
  DerivedAnnotation(AbstractDerivedAnnotation<?> annotation) {
    mAnnotation = annotation;
    mType = getCompatibleType(annotation.getField().getType());
  }

  /**
   * @param annotationName name of the annotation
   */
  public DerivedAnnotation(String annotationName) {
    this(DerivedAnnotations.valueOf(annotationName).getAnnotation());
  }

  protected static AnnotationDataType getCompatibleType(MetaType mt) {
    switch (mt) {
      case INTEGER:
        return AnnotationDataType.INTEGER;
      case FLOAT:
        return AnnotationDataType.DOUBLE;
      case STRING:
      case CHARACTER:
        return AnnotationDataType.STRING;
      case FLAG:
        return AnnotationDataType.BOOLEAN;
      default:
        return null;
    }
  }

  @Override
  public String getName() {
    return "DERIVED-" + mAnnotation.getName();
  }

  @Override
  public AnnotationDataType getType() {
    return mType;
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
