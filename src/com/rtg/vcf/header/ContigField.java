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

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Pattern;

import com.rtg.util.StringUtils;
import com.rtg.util.Utils;

/**
 *
 */
public class ContigField implements IdField<ContigField> {

  private static final Pattern CONTIG_LINE_PATTERN = Pattern.compile("^##contig=<(.+)>$");

  private final String mId;
  private final Integer mLength;
  private final Map<String, String> mValues;

  /**
   * @param line filter line from <code>VCF</code> file
   */
  public ContigField(String line) {
    final LinkedHashMap<String, String> temp = VcfHeader.parseMetaLine(line, CONTIG_LINE_PATTERN);
    mId = temp.get("ID");
    if (mId == null) {
      throw new IllegalArgumentException("Expected ID field on line: " + line);
    }
    temp.remove("ID");
    if (temp.containsKey("length")) {
      mLength = Integer.valueOf(temp.get("length"));
      temp.remove("length");
    } else {
      mLength = null;
    }
    mValues = temp;
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof ContigField)) {
      return false;
    }
    final ContigField other = (ContigField) obj;
    return mId.equals(other.mId) && mValues.equals(other.mValues)
      && ((mLength == null && other.mLength == null)
          || (mLength != null && mLength.equals(other.mLength)));
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(Utils.pairHash(mId.hashCode(), mValues.hashCode()), mLength != null ? mLength : 0);
  }

  @Override
  public ContigField superSet(ContigField other) {
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
  public Integer getLength() {
    return mLength;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append(VcfHeader.CONTIG_STRING).append("=<ID=").append(mId);
    if (mLength != null) {
      sb.append(",length=").append(mLength);
    }
    for (final Map.Entry<String, String> entry : mValues.entrySet()) {
      sb.append(",").append(entry.getKey()).append("=");
      sb.append(StringUtils.smartQuote(entry.getValue()));
    }
    sb.append(">");
    return sb.toString();
  }
}
