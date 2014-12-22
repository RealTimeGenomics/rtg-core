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
 * Class to encapsulate information of an alt meta information line in <code>VCF</code>
 */
public class AltField implements IdField<AltField> {

  private static final Pattern ALT_LINE_PATTERN = Pattern.compile("^##ALT=<(.+)>$");

  private final String mId;
  private final String mDescription;

  /**
   * @param line filter line from <code>VCF</code> file
   */
  public AltField(String line) {
    final HashMap<String, String> temp = VcfHeader.parseMetaLine(line, ALT_LINE_PATTERN);
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
  public AltField(String id, String description) {
    mId = id;
    mDescription = description;
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof AltField)) {
      return false;
    }
    final AltField other = (AltField) obj;
    return mId.equals(other.mId) && mDescription.equals(other.mDescription);
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mId.hashCode(), mDescription.hashCode());
  }

  @Override
  public AltField superSet(AltField other) {
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
    return VcfHeader.ALT_STRING + "=<ID=" + mId + ",Description=\"" + mDescription + "\">";
  }

}
