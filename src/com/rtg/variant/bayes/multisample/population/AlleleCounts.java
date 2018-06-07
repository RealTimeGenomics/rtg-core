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

package com.rtg.variant.bayes.multisample.population;

import java.util.Map;
import java.util.Set;

/**
 * Allele counts for population variants.
 */
public class AlleleCounts {
  private final int mPosition;
  private final Map<String, Integer> mCounts;
  private final String mReferenceAllele;
  private final boolean mIsComplex;
  private final int mMaxLength;

  /**
   *
   * @param position zero based position of alleles
   * @param counts map of allele to counts
   * @param referenceAllele the reference allele
   */
  public AlleleCounts(int position, Map<String, Integer> counts, String referenceAllele) {
    if (counts == null) {
      throw new NullPointerException();
    }
    mPosition = position;
    mCounts = counts;
    mReferenceAllele = referenceAllele;
    final int rlength = referenceAllele.length();
    boolean complex = rlength != 1;
    int max = rlength;
    for (final String s : counts.keySet()) {
      final int length = s.length();
      if (length > max) {
        max = length;
      }
      complex |= length != 1;
    }
    mIsComplex = complex;
    mMaxLength = max;
  }

  /**
   * @return zero based template position of variants
   */
  public int position() {
    return mPosition;
  }

  /**
   * Return count for a given allele.
   * @param allele an allele
   * @return allele count, or -1 if allele does not exist
   */
  public int count(String allele) {
    return mCounts.getOrDefault(allele, -1);
  }

  /**
   * Returns the set of alleles that there are counts for.
   * @return set of alleles
   */
  public Set<String> allelesSeen() {
    return mCounts.keySet();
  }

  @Override
  public String toString() {
    return "(" + mPosition + " " + mCounts + ")";
  }

  /**
   * @return true if one or more alleles represents a complex call (not a SNP)
   */
  public boolean isComplex() {
    return mIsComplex;
  }

  /**
   * @return the maximum length of the reference and all the alleles.
   */
  public int maxLength() {
    return mMaxLength;
  }

  public String getReferenceAllele() {
    return mReferenceAllele;
  }

  /**
   * Returns the length of reference that the set of alleles apply to
   * @return number of bases in reference
   */
  public int refLength() {
    return mReferenceAllele.length();
  }
}
