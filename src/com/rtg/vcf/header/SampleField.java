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
import java.util.Locale;
import java.util.regex.Pattern;

import com.rtg.reference.Sex;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;

/**
 * Class to encapsulate a sample line in <code>VCF</code>
 */
public class SampleField implements IdField<SampleField> {

  private static final Pattern SAMPLE_LINE_PATTERN = Pattern.compile("^##SAMPLE=<(.+)>$");

  private final String mId;
  private final String mGenomes;
  private final String mMixture;
  private final String mDescription;
  private final Sex mSex;

  /**
   * @param line sample line
   */
  public SampleField(String line) {
    final HashMap<String, String> temp = VcfHeader.parseMetaLine(line, SAMPLE_LINE_PATTERN);
    VcfHeader.checkRequiredMetaKeys(temp, line, "ID");
    mId = temp.get("ID");
    mGenomes = temp.get("Genomes");
    mMixture = temp.get("Mixture");
    mDescription = temp.get("Description");
    final String sexStr = temp.get("Sex");
    if (sexStr != null) {
      mSex = Sex.valueOf(sexStr.toUpperCase(Locale.getDefault()));
    } else {
      mSex = null;
    }
  }

  /**
   * @return the ID field
   */
  @Override
  public String getId() {
    return mId;
  }

  /**
   * @return genomes field
   */
  public String getGenomes() {
    return mGenomes;
  }

  /**
   * @return mixture field
   */
  public String getMixture() {
    return mMixture;
  }

  /**
   * @return description field
   */
  public String getDescription() {
    return mDescription;
  }

  public Sex getSex() {
    return mSex;
  }

  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof SampleField)) {
      return false;
    }
    final SampleField other = (SampleField) obj;
    return mId.equals(other.mId) && Utils.equals(mGenomes, other.mGenomes)
      && Utils.equals(mMixture, other.mMixture) && Utils.equals(mDescription, other.mDescription);
  }

  @Override
  public int hashCode() {
    final int gHash = mGenomes != null ? mGenomes.hashCode() : 0;
    final int mxHash = mMixture != null ? mMixture.hashCode() : 0;
    final int dHash = mDescription != null ? mDescription.hashCode() : 0;
    final int sxHash = mSex != null ? mSex.hashCode() : 0;
    return Utils.pairHashContinuous(mId.hashCode(), gHash, mxHash, dHash, sxHash);
  }

  @Override
  public SampleField superSet(SampleField other) {
    return equals(other) ? this : null;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append(VcfHeader.SAMPLE_STRING).append("=<ID=").append(mId);
    if (mGenomes != null) {
      sb.append(",Genomes=").append(mGenomes);
    }
    if (mMixture != null) {
      sb.append(",Mixture=").append(mMixture);
    }
    if (mDescription != null) {
      sb.append(",Description=").append(StringUtils.dumbQuote(mDescription));
    }
    if (mSex != null && (mSex == Sex.MALE || mSex == Sex.FEMALE)) {
      sb.append(",Sex=").append(mSex.toString());
    }
    sb.append(">");
    return sb.toString();
  }
}
