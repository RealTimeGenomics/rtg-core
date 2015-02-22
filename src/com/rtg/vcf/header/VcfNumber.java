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

import com.rtg.util.Utils;

/**
 * Encapsulate value of number field from <code>VCF</code> meta information lines
 */
public class VcfNumber {
  private final int mNumber;
  private final VcfNumberType mType;

  /** Singleton for variable numbered variant types */
  public static final VcfNumber DOT = new VcfNumber(".");
  /** Singleton for alt allele valued variant types */
  public static final VcfNumber ALTS = new VcfNumber("A");
  /** Singleton for all allele valued variant types */
  public static final VcfNumber REF_ALTS = new VcfNumber("R");
  /** Singleton for genotype valued variant types */
  public static final VcfNumber GENOTYPES = new VcfNumber("G");

  /**
   * @param number number field from a meta line
   */
  public VcfNumber(String number) {
    switch (number) {
      case "A":
        mType = VcfNumberType.ALTS;
        mNumber = -1;
        break;
      case "R":
        mType = VcfNumberType.REF_ALTS;
        mNumber = -1;
        break;
      case "G":
        mType = VcfNumberType.GENOTYPES;
        mNumber = -1;
        break;
      case ".":
        mType = VcfNumberType.UNKNOWN;
        mNumber = -1;
        break;
      default:
        mType = VcfNumberType.INTEGER;
        mNumber = Integer.parseInt(number);
        break;
    }
  }
  @Override
  public boolean equals(Object obj) {
    if (!(obj instanceof VcfNumber)) {
      return false;
    }
    final VcfNumber other = (VcfNumber) obj;
    switch (mType) {
      case ALTS:
      case GENOTYPES:
      case REF_ALTS:
      case UNKNOWN:
        return mType == other.mType;
      case INTEGER:
        return other.mType == VcfNumberType.INTEGER && mNumber == other.mNumber;
      default:
        throw new IllegalArgumentException("Type: " + mType + " is not supported");
    }
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mType.ordinal(), mNumber);
  }

  /**
   * @return the numbers type
   */
  public VcfNumberType getNumberType() {
    return mType;
  }

  /**
   * @return if type is {@link VcfNumberType#INTEGER} then its value, otherwise -1
   */
  public int getNumber() {
    return mNumber;
  }

  @Override
  public String toString() {
    return mType.isFixedString() ? mType.toString() : Integer.toString(mNumber);
  }

}
