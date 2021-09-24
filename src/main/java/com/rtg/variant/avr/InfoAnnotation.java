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
