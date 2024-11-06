/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
   * @return reference allele frequency, or 1.0 if all counts are zero
   */
  public double refAlleleFrequency() {
    int sum = 0;
    int refcount = 0;
    for (Map.Entry<String, Integer> entry : mCounts.entrySet()) {
      final int c = entry.getValue();
      sum += c;
      if (mReferenceAllele.equals(entry.getKey())) {
        refcount = c;
      }
    }
    return sum == 0 ? 1.0 : (double) refcount / sum;
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
