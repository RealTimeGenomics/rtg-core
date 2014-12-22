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
package com.rtg.vcf.header;

import java.util.HashMap;
import java.util.regex.Pattern;

import com.rtg.util.Utils;

/**
 * Class to encapsulate information of a info meta information line in <code>VCF</code>
 */
public class InfoField implements IdField<InfoField> {
  private static final Pattern INFO_LINE_PATTERN = Pattern.compile("^##INFO=<(.+)>$");

  private final String mId;
  private final MetaType mType;
  private final VcfNumber mNumber;
  private final String mDescription;

  /**
   * @param line info line from <code>VCF</code> file
   */
  public InfoField(String line) {
    final HashMap<String, String> temp = VcfHeader.parseMetaLine(line, INFO_LINE_PATTERN);
    mId = temp.get("ID");
    if (mId == null) {
      throw new IllegalArgumentException("Expected ID field on line: " + line);
    }
    mType = MetaType.parseValue(temp.get("Type"));
    mNumber = new VcfNumber(temp.get("Number"));
    mDescription = temp.get("Description");
  }

  /**
   * Makes a VCF InfoField
   * @param id the info field identifier
   * @param type the type of value
   * @param number the specifier for the number of occurrences
   * @param description the field description
   */
  public InfoField(String id, MetaType type, VcfNumber number, String description) {
    mId = id;
    mType = type;
    mNumber = number;
    mDescription = description;
  }

  @Override
  public boolean equals(Object obj) {
    return mostlyEquals(obj) && mType == ((InfoField) obj).mType;
  }

  private boolean mostlyEquals(Object obj) {
    if (!(obj instanceof InfoField)) {
      return false;
    }
    final InfoField other = (InfoField) obj;
    return mId.equals(other.mId) && mNumber.equals(other.mNumber) && mDescription.equals(other.mDescription);
  }

  @Override
  public InfoField superSet(InfoField other) {
    if (!mostlyEquals(other)) {
      return null;
    }
    if (mType.isSuperSet(other.mType)) {
      return this;
    } else if (other.mType.isSuperSet(mType)) {
      return other;
    }
    return null;
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(Utils.pairHash(Utils.pairHash(mId.hashCode(), mNumber.hashCode()), mType.ordinal()), mDescription.hashCode());
  }

  /**
   * @return ID field
   */
  @Override
  public String getId() {
    return mId;
  }

  /**
   * @return number field
   */
  public VcfNumber getNumber() {
    return mNumber;
  }

  /**
   * @return type field
   */
  public MetaType getType() {
    return mType;
  }

  /**
   * @return description field
   */
  public String getDescription() {
    return mDescription;
  }

  @Override
  public String toString() {
    return VcfHeader.INFO_STRING + "=<ID=" + mId + ",Number=" + mNumber.toString() + ",Type=" + mType.toString() + ",Description=\"" + mDescription + "\">";
  }

}
