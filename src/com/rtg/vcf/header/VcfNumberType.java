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

/**
 * Enum for type of number in <code>VCF</code> meta lines
 */
public enum VcfNumberType {

  /** if specified number of values applies to field */
  INTEGER(false, ""),
  /** if one value per alternate allele */
  ALTS(true, "A"),
  /** if one value per all possible alleles (ref and alts) */
  REF_ALTS(true, "R"),
  /** if one value per genotype */
  GENOTYPES(true, "G"),
  /** if number of values varies or is unknown */
  UNKNOWN(true, ".");

  private final boolean mFixed;
  private final String mToString;

  private VcfNumberType(boolean fixedString, String str) {
    mFixed = fixedString;
    mToString = str;
  }

  /**
   * @return true if {@link VcfNumberType#toString()} will return numbers output value, false if should use {@link VcfNumber#toString()} instead
   */
  public boolean isFixedString() {
    return mFixed;
  }

  @Override
  public String toString() {
    return mFixed ? mToString : super.toString();
  }

}
