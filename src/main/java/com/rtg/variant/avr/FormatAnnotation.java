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
import java.util.List;

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * An annotation extracted from an existing FORMAT field in a VCF file.
 */
public class FormatAnnotation implements Annotation {
  private final String mFieldName;
  private final AnnotationDataType mType;

  /**
   * Creates an FORMAT annotation for the given name and type.
   * Annotation types must have a <code>valueOf(String)</code> method, or be a <code>String</code>.
   * @param fieldName annotation field name
   * @param type annotation type
   */
  public FormatAnnotation(String fieldName, AnnotationDataType type) {
    if (fieldName == null) {
      throw new NullPointerException("null name given.");
    }
    if (type == null) {
      throw new NullPointerException("null type given.");
    }
    mFieldName = fieldName;
    if (type != AnnotationDataType.INTEGER
        && type != AnnotationDataType.DOUBLE
        && type != AnnotationDataType.STRING) {
      throw new IllegalArgumentException("Invalid format type: " + type);
    }
    mType = type;
  }

  @Override
  public String getName() {
    return "FORMAT-" + mFieldName;
  }

  @Override
  public AnnotationDataType getType() {
    return mType;
  }

  @Override
  public Object getValue(VcfRecord record, int sampleNumber) {
    final List<String> sampleValues = record.getFormat(mFieldName);
    final Object value;
    if (sampleValues == null || sampleValues.size() < (sampleNumber + 1)) {
      value = null;
    } else {
      final String sValue = sampleValues.get(sampleNumber);
      if (VcfRecord.MISSING.equals(sValue)) {
        value = null;
      } else {
        if (sValue.indexOf(',') >= 0) {
          throw new IllegalArgumentException("We don't support multi-value FORMAT fields: " + mFieldName);
        }
        value = mType.stringToObjectOfType(sValue);
      }
    }
    return value;
  }

  @Override
  public String checkHeader(VcfHeader header) {
    for (FormatField field : header.getFormatLines()) {
      if (field.getId().equals(mFieldName)) {
        if (!getType().isMetaTypeCompatible(field.getType())) {
          return "FORMAT field type mismatch: " + mFieldName + " : " + field.getType() + " != " + getType();
        }
        final VcfNumber number = field.getNumber();
        switch (number.getNumberType()) {
          case INTEGER:
            if (number.getNumber() != 1) {
              return "FORMAT field too many values: " + mFieldName + " : " + number.getNumber() + " != 1";
            }
            break;
          default:
            return "FORMAT field arbitary number of values: " + mFieldName;
        }
        return null; // all tests for this field passed
      }
    }
    return "FORMAT field does not exist: " + mFieldName;
  }

  @Override
  public void save(DataOutputStream dos) throws IOException {
    dos.writeInt(AnnotationLoader.AnnotationType.FORMAT.ordinal());
    dos.writeUTF(mFieldName);
    dos.writeInt(getType().ordinal());
  }

}
