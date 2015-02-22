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
 * Class to encapsulate information of a format meta information line in <code>VCF</code>
 */
public class FormatField implements IdField<FormatField> {

  private static final Pattern FORMAT_LINE_PATTERN = Pattern.compile("^##FORMAT=<(.+)>$");

  private final String mId;
  private final MetaType mType;
  private final VcfNumber mNumber;
  private final String mDescription;

  /**
   * @param line format line
   */
  public FormatField(String line) {
    final HashMap<String, String> temp = VcfHeader.parseMetaLine(line, FORMAT_LINE_PATTERN);
    VcfHeader.checkRequiredMetaKeys(temp, line, "ID", "Type", "Number", "Description");
    mId = temp.get("ID");
    mType = MetaType.parseValue(temp.get("Type"));
    mNumber = new VcfNumber(temp.get("Number"));
    mDescription = temp.get("Description");
  }

  /**
   * Makes a VCF FormatField
   * @param id the format field identifier
   * @param type the type of value
   * @param number the specifier for the number of occurrences
   * @param description the field description
   */
  public FormatField(String id, MetaType type, VcfNumber number, String description) {
    mId = id;
    mType = type;
    mNumber = number;
    mDescription = description;
  }


  /**
   * @return the ID field
   */
  @Override
  public String getId() {
    return mId;
  }

  @Override
  public boolean equals(Object obj) {
    return mostlyEquals(obj) && mType == ((FormatField) obj).mType;
  }

  private boolean mostlyEquals(Object obj) {
    if (!(obj instanceof FormatField)) {
      return false;
    }
    final FormatField other = (FormatField) obj;
    return mId.equals(other.mId) && mNumber.equals(other.mNumber) && mDescription.equals(other.mDescription);
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(Utils.pairHash(Utils.pairHash(mId.hashCode(), mNumber.hashCode()), mType.ordinal()), mDescription.hashCode());
  }

  @Override
  public FormatField superSet(FormatField other) {
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
    return VcfHeader.FORMAT_STRING + "=<ID=" + mId + ",Number=" + mNumber.toString() + ",Type=" + mType.toString() + ",Description=\"" + mDescription + "\">";
  }

}
