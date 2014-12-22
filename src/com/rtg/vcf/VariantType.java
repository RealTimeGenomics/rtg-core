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
package com.rtg.vcf;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Enumeration of Variant types
 */
@TestClass("com.rtg.variant.VariantStatisticsTest")
public enum VariantType {
  //Ordered least to most precedence

  /** Equal to reference */
  UNCHANGED(false, false),
  /** Single nucleotide polymorphism */
  SNP(false, false),
  /** Multiple nucleotide polymorphism */
  MNP(false, false),
  /** Deletion with respect to reference */
  DELETION(true, false),
  /** Insertion with respect to reference */
  INSERTION(true, false),
  /** Complex indel with deletions and insertions */
  INDEL(true, false),
  /** Structural variant break end */
  SV_BREAKEND(false, true),
  /** Structural variant symbolic value */
  SV_SYMBOLIC(false, true);

  private final boolean mIsIndelType;
  private final boolean mIsSvType;

  VariantType(boolean isIndelType, boolean isSvType) {
    mIsIndelType = isIndelType;
    mIsSvType = isSvType;
  }

  public boolean isIndelType() {
    return mIsIndelType;
  }

  public boolean isSvType() {
    return mIsSvType;
  }

  /**
   * If the allele is a symbolic structural variant, return the appropriate VariantType
   * @param allele the alt allele from VCF.
   * @return <code>VariantType.SV_SYMBOLIC</code> or <code>VariantType.SV_BREAKEND</code>, or null if not a symbolic variant
   */
  public static VariantType getSymbolicAlleleType(String allele) {
    if (allele.length() > 0) {
      if (allele.charAt(0) == '<') {
        return VariantType.SV_SYMBOLIC;
      } else {
        for (int i = 0; i < allele.length(); i++) {
          final char c = allele.charAt(i);
          if ((c == '[') || (c == ']')) {
            return VariantType.SV_BREAKEND;
          }
        }
      }
    }
    return null;
  }

  /**
   * Determine the type of variant given a reference allele and an alternative allele
   * @param refAllele the reference allele
   * @param altAllele the alternative allele
   * @return the type of variant
   */
  public static VariantType getType(String refAllele, String altAllele) {
    if (refAllele.equals(altAllele)) {
      return VariantType.UNCHANGED;
    } else if (refAllele.length() == 1 && altAllele.length() == 1) {
      return VariantType.SNP;
    }
    final VariantType svType = getSymbolicAlleleType(altAllele);
    if (svType != null) {
      return svType;
    }
    if (refAllele.length() == altAllele.length()) {
      return VariantType.MNP;
    } else if (isInsertionOrDeletion(refAllele, altAllele)) {
      if (refAllele.length() < altAllele.length()) {
        return VariantType.INSERTION;
      } else {
        return VariantType.DELETION;
      }
    } else {
      return VariantType.INDEL;
    }
  }

  private static boolean isInsertionOrDeletion(String ref, String pred) {
    //Only call if have already ruled out Unchanged, SNP and MNP
    if (ref.length() == 0 || pred.length() == 0) {
      return true;
    }
    final boolean isInsert = ref.length() < pred.length();
    final String seq1 = isInsert ? ref : pred;
    final String seq2 = isInsert ? pred : ref;
    final int seq1Len = seq1.length();
    final int seq2Len = seq2.length();
    int left;
    for (left = 0; left < seq1Len; left++) {
      if (seq1.charAt(left) != seq2.charAt(left)) {
        break;
      }
    }
    if (left == seq1Len) {
      return true;
    }
    int right;
    for (right = 0; right < seq1Len; right++) {
      if (seq1.charAt(seq1Len - right - 1) != seq2.charAt(seq2Len - right - 1)) {
        break;
      }
    }
    return left + right >= seq1Len;
  }
}
