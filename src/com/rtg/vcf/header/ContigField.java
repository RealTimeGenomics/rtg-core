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
    VcfHeader.checkRequiredMetaKeys(temp, line, "ID");
    mId = temp.get("ID");
    temp.remove("ID");
    if (temp.containsKey("length")) {
      mLength = Integer.valueOf(temp.get("length"));
      temp.remove("length");
    } else {
      mLength = null;
    }
    mValues = temp;
  }

  ContigField(String id, Integer length) {
    mId = id;
    mLength = length;
    mValues = new LinkedHashMap<>();
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof ContigField)) {
      return false;
    }
    final ContigField other = (ContigField) obj;
    return mostlyEquals(other)
      && ((mLength == null) == (other.mLength == null))
      && (mLength == null || mLength.equals(other.mLength))
      && mValues.equals(other.mValues);
  }

  // True if no conflicting fields preventing merge
  private boolean mostlyEquals(ContigField other) {
    if (!mId.equals(other.mId)) {
      return false;
    }
    if (mLength != null && other.mLength != null && !mLength.equals(other.mLength)) {
      return false;
    }
    for (Map.Entry<String, String> entry : mValues.entrySet()) {
      final String val = other.mValues.get(entry.getKey());
      if (val != null && !val.equals(entry.getValue())) {
        return false;
      }
    }
    return true;
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(Utils.pairHash(mId.hashCode(), mValues.hashCode()), mLength != null ? mLength : 0);
  }

  @Override
  public ContigField superSet(ContigField other) {
    if (!mostlyEquals(other)) {
      return null;
    }
    if (equals(other)) {
      return this;
    }
    final ContigField result = new ContigField(mId, mLength != null ? mLength : other.mLength);
    result.mValues.putAll(mValues);
    result.mValues.putAll(other.mValues);
    return result;
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
