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
import java.util.ArrayList;

import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.annotation.AnnotationDataType;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;

/**
 * Annotations for FORMAT fields in VCF files.
 *
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
    final ArrayList<String> sampleValues = record.getFormatAndSample().get(mFieldName);
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
          return "FORMAT field type mismatch: " + mFieldName + " : " + field.getType() + " != " + getType().toString();
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
