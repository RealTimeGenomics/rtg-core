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
package com.rtg.variant.util;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;

import com.rtg.vcf.VcfUtils;

/**
 * Represents a called genotype, with numeric allele coding according to VCF (missing is encoded as -1). Alleles are in sorted order.
 */
class Genotype {

  static final GenotypeComparator GENOTYPE_COMPARATOR = new GenotypeComparator();

  private final int[] mAlleles;

  Genotype(String gt) {
    this(VcfUtils.splitGt(gt));
  }

  Genotype(int[] gt) {
    mAlleles = gt;
    Arrays.sort(mAlleles);
  }

  public int length() {
    return mAlleles.length;
  }

  public int get(int index) {
    return mAlleles[index];
  }

  public boolean contains(int allele) {
    for (int a : mAlleles) {
      if (a == allele) {
        return true;
      }
    }
    return false;
  }

  public boolean homozygous() {
    return mAlleles.length > 0 && mAlleles[0] == mAlleles[mAlleles.length - 1];
  }

  public String toString() {
    final StringBuilder sb = new StringBuilder();
    boolean first = true;
    for (int j : mAlleles) {
      if (first) {
        sb.append(Integer.toString(j));
        first = false;
      } else {
        sb.append("/").append(Integer.toString(j));
      }
    }
    return sb.toString();
  }

  public boolean equals(Object o) {
    if (o == null || !(o instanceof Genotype)) {
      return false;
    }
    return GENOTYPE_COMPARATOR.compare(this, (Genotype) o) == 0;
  }

  public int hashCode() {
    int res = 0;
    for (int h : mAlleles) {
      res += h;
    }
    return res;
  }

  // If this modified to be non-thread safe, then make sure to remove the static instance
  static class GenotypeComparator implements Comparator<Genotype>, Serializable {
    @Override
    public int compare(Genotype a, Genotype b) {
      for (int i = 0; i < a.mAlleles.length && i < b.mAlleles.length; i++) {
        if (a.mAlleles[i] < b.mAlleles[i]) {
          return -1;
        } else if (a.mAlleles[i] > b.mAlleles[i]) {
          return 1;
        }
      }
      if (a.mAlleles.length < b.mAlleles.length) {
        return -1;
      } else if (a.mAlleles.length > b.mAlleles.length) {
        return 1;
      }
      return 0;
    }
  }
}
