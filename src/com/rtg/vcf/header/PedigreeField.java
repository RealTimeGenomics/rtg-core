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
import java.util.Map;
import java.util.regex.Pattern;

import com.rtg.util.Utils;

/**
 * Class to encapsulate a pedigree line in <code>VCF</code>
 */
public class PedigreeField {

  //Family pedigree relationship fields
  private static final String PEDIGREE_FATHER = "Father";
  private static final String PEDIGREE_MOTHER = "Mother";
  private static final String PEDIGREE_CHILD = "Child";

  //Clonal pedigree relationship fields
  private static final String PEDIGREE_DERIVED = "Derived";
  private static final String PEDIGREE_ORIGINAL = "Original";

  private static final Pattern PEDIGREE_LINE_PATTERN = Pattern.compile("^##PEDIGREE=<(.+)>$");
  private final HashMap<String, String> mSamples;

  /**
   * @param line pedigree line
   */
  public PedigreeField(String line) {
    mSamples = VcfHeader.parseMetaLine(line, PEDIGREE_LINE_PATTERN);
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof PedigreeField)) {
      return false;
    }
    final PedigreeField other = (PedigreeField) obj;
    if (mSamples.size() != other.mSamples.size()) {
        return false;
    }
    for (final String key : mSamples.keySet()) {
      if (!mSamples.get(key).equals(other.mSamples.get(key))) {
        return false;
      }
    }
    return true;
  }

  @Override
  public int hashCode() {
    int h = 0;
    for (final Map.Entry<String, String> e : mSamples.entrySet()) {
      h = Utils.pairHash(h, Utils.pairHash(e.getKey().hashCode(), e.getValue().hashCode()));
    }
    return h;
  }

  /**
   * @return the child sample
   */
  public String getChild() {
    return mSamples.get(PEDIGREE_CHILD);
  }

  /**
   * @return the mother sample
   */
  public String getMother() {
    return mSamples.get(PEDIGREE_MOTHER);
  }

  /**
   * @return the father sample
   */
  public String getFather() {
    return mSamples.get(PEDIGREE_FATHER);
  }

  /**
   * @return the derived sample
   */
  public String getDerived() {
    return mSamples.get(PEDIGREE_DERIVED);
  }

  /**
   * @return the original original sample
   */
  public String getOriginal() {
    return mSamples.get(PEDIGREE_ORIGINAL);
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder(VcfHeader.PEDIGREE_STRING).append("=<");
    boolean first = true;
    for (final Map.Entry<String, String> e : mSamples.entrySet()) {
      if (first) {
        first = false;
      } else {
        sb.append(",");
      }
      sb.append(e.getKey()).append("=").append(e.getValue());
    }
    sb.append(">");
    return sb.toString();
  }
}
