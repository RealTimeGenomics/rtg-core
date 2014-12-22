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
 * Class to encapsulate information of a filter meta information line in <code>VCF</code>
 */
public class FilterField implements IdField<FilterField> {

  private static final Pattern FILTER_LINE_PATTERN = Pattern.compile("^##FILTER=<(.+)>$");

  private final String mId;
  private final String mDescription;

  /**
   * @param line filter line from <code>VCF</code> file
   */
  public FilterField(String line) {
    final HashMap<String, String> temp = VcfHeader.parseMetaLine(line, FILTER_LINE_PATTERN);
    mId = temp.get("ID");
    if (mId == null) {
      throw new IllegalArgumentException("Expected ID field on line: " + line);
    }
    mDescription = temp.get("Description");
  }

  /**
   * Makes a VCF FilterField
   * @param id the filter field identifier
   * @param description the field description
   */
  public FilterField(String id, String description) {
    mId = id;
    mDescription = description;
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof FilterField)) {
      return false;
    }
    final FilterField other = (FilterField) obj;
    return mId.equals(other.mId) && mDescription.equals(other.mDescription);
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mId.hashCode(), mDescription.hashCode());
  }

  @Override
  public FilterField superSet(FilterField other) {
    return equals(other) ? this : null;
  }

  /**
   * @return the ID field
   */
  @Override
  public String getId() {
    return mId;
  }

  /**
   * @return the description field
   */
  public String getDescription() {
    return mDescription;
  }

  @Override
  public String toString() {
    return VcfHeader.FILTER_STRING + "=<ID=" + mId + ",Description=\"" + mDescription + "\">";
  }

}
