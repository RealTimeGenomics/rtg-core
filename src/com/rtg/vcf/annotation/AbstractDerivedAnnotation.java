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

import java.util.HashSet;
import java.util.Set;

import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.VcfHeader;

/**
 * Abstract class to use when implementing a derived annotation.
 */
public abstract class AbstractDerivedAnnotation implements VcfAnnotation {

  private final String mName;
  private final String mDescription;
  private final AnnotationDataType mType;

  /**
   * @param name the attribute name
   * @param description the attribute description
   * @param type the attribute type
   */
  public AbstractDerivedAnnotation(String name, String description, AnnotationDataType type) {
    mName = name;
    mDescription = description;
    mType = type;
  }

  @Override
  public String getName() {
    return mName;
  }

  @Override
  public AnnotationDataType getType() {
    return mType;
  }

  /**
   * Get the VCF header description for the annotation.
   * @return the description of the annotation.
   */
  public String getDescription() {
    return mDescription;
  }

  protected String checkHeader(VcfHeader header, String[] infoFields, String[] formatFields) {
    final Set<String> infoHeaderIds = new HashSet<>();
    final Set<String> formatHeaderIds = new HashSet<>();
    if (header != null) {
      for (final InfoField field : header.getInfoLines()) {
        infoHeaderIds.add(field.getId());
      }
      for (final FormatField field : header.getFormatLines()) {
        formatHeaderIds.add(field.getId());
      }
    }
    final StringBuilder missingInfos = new StringBuilder();
    if (infoFields != null) {
      for (final String info : infoFields) {
        if (!infoHeaderIds.contains(info)) {
          missingInfos.append(' ').append(info);
        }
      }
    }
    final StringBuilder missingFormats = new StringBuilder();
    if (formatFields != null) {
      for (final String format : formatFields) {
        if (!formatHeaderIds.contains(format)) {
          missingFormats.append(' ').append(format);
        }
      }
    }
    final StringBuilder sb = new StringBuilder();
    if (missingInfos.length() > 0 || missingFormats.length() > 0) {
      sb.append("Derived annotation ").append(mName).append(" missing required fields in VCF header");
      if (missingInfos.length() > 0) {
        sb.append(" (INFO fields:").append(missingInfos.toString()).append(')');
      }
      if (missingFormats.length() > 0) {
        sb.append(" (FORMAT fields:").append(missingFormats.toString()).append(')');
      }
    }
    return sb.length() > 0 ? sb.toString() : null;
  }

}
