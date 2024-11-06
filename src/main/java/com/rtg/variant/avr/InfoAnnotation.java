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
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * An annotation extracted from an existing INFO field in a VCF file.
 */
public class InfoAnnotation implements Annotation {

  private final String mFieldName;
  private final AnnotationDataType mType;

  /**
   * Creates an INFO annotation for the given name and type.
   * Annotation types must have a <code>valueOf(String)</code> method, or be a <code>String</code>.
   * If null if given for the type, it is assumed to be a flag (true if present, false otherwise).
   * @param fieldName INFO field name
   * @param type annotation type
   */
  public InfoAnnotation(String fieldName, AnnotationDataType type) {
    if (fieldName == null) {
      throw new NullPointerException("null name given.");
    }
    mFieldName = fieldName;
    if (type == null) {
      throw new NullPointerException("null type given");
    }
    if (type != AnnotationDataType.BOOLEAN
        && type != AnnotationDataType.INTEGER
        && type != AnnotationDataType.DOUBLE
        && type != AnnotationDataType.STRING) {
      throw new IllegalArgumentException("Invalid info type: " + type);
    }
    mType = type;
  }

  @Override
  public String getName() {
    return "INFO-" + mFieldName;
  }

  @Override
  public AnnotationDataType getType() {
    return mType;
  }

  /**
   * Returns the value for an annotation from a given {@link VcfRecord}.
   * For sample specific annotations, the value in the first sample is returned.
   * @param record record to extract value from.
   * @param sampleNumber not used for INFO annotations
   * @return value for this annotation from the given record.
   * @throws IllegalArgumentException cannot convert VCF record value to annotation type
   */
  @Override
  public Object getValue(VcfRecord record, int sampleNumber) {
    final String[] infoValues = record.getInfoSplit(mFieldName);
    final Object value;
    if (infoValues == null) {
      value = getType() == AnnotationDataType.BOOLEAN ? Boolean.FALSE : null;
    } else if (infoValues.length == 0) {
      if (getType() != AnnotationDataType.BOOLEAN) {
        throw new IllegalArgumentException("Value expected for type: " + getType());
      }
      value = Boolean.TRUE;
    } else if (infoValues.length > 1) {
      throw new IllegalArgumentException("We don't support multi-value INFO fields: " + mFieldName);
    } else {
      value = mType.stringToObjectOfType(infoValues[0]);
    }
    return value;
  }

  @Override
  public String checkHeader(VcfHeader header) {
    for (InfoField field : header.getInfoLines()) {
      if (field.getId().equals(mFieldName)) {
        if (!getType().isMetaTypeCompatible(field.getType())) {
          return "INFO field type mismatch: " + mFieldName + " : " + field.getType() + " != " + getType();
        }
        final VcfNumber number = field.getNumber();
        switch (number.getNumberType()) {
          case INTEGER:
            if (number.getNumber() > 1) {
              return "INFO field too many values: " + mFieldName + " : " + number.getNumber() + " > 1";
            }
            break;
          default:
            return "INFO field arbitary number of values: " + mFieldName;
        }
        return null; // all tests for this field passed
      }
    }
    return "INFO field does not exist: " + mFieldName;
  }

  @Override
  public void save(DataOutputStream dos) throws IOException {
    dos.writeInt(AnnotationLoader.AnnotationType.INFO.ordinal());
    dos.writeUTF(mFieldName);
    dos.writeInt(getType().ordinal());
  }

}
