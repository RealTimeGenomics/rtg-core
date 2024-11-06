/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
