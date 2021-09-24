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

package com.rtg.segregation;

import com.rtg.reference.Ploidy;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;

/**
 * Single diploid genotype that only knows about two allele numbers from VCF.
 */
class GType extends IntegralAbstract {

  private final int mA;
  private final int mB;
  private final Ploidy mPloidy;

  GType(int a, int b, Ploidy ploidy) {
    mA = a;
    mB = b;
    mPloidy = ploidy;
  }

  GType(final String str, final Ploidy ploidy) throws MismatchingPloidyException {
    assert ploidy != null;
    mPloidy = ploidy;
    if (ploidy == Ploidy.NONE) {
      mA = -1;
      mB = -1;
    } else {
      if (VcfRecord.MISSING.equals(str)) {
        mA = 0;
        if (ploidy == Ploidy.DIPLOID) {
          mB = 0;
        } else {
          mB = -1;
        }
      } else {
        final int[] split = VcfUtils.splitGt(str);
        final int a = split[0];
        if (split.length > 1) {
          if (ploidy != Ploidy.DIPLOID) {
            throw new MismatchingPloidyException(str, ploidy);
          }
          final int b = split[1];
          mA = Math.min(a, b);
          mB = Math.max(a, b);
        } else {
          if (ploidy != Ploidy.HAPLOID && ploidy != Ploidy.POLYPLOID) {
            throw new MismatchingPloidyException(str, ploidy);
          }
          mA = a;
          mB = -1;
        }
      }
    }
  }

  int a() {
    return mA;
  }

  int b() {
    return mB;
  }

  Ploidy ploidy() {
    return mPloidy;
  }

  /**
   * @return true iff only single allele.
   */
  boolean isSingleAllele() {
    return mA == mB || mB == -1;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(0 <= mA && mA <= mB);
    return true;
  }

  @Override
  public void toString(StringBuilder sb) {
    if (mPloidy == Ploidy.NONE) {
      sb.append("NONE");
    } else {
      sb.append(mA);
      if (mB > -1) {
        sb.append("_").append(mB);
      }
    }
  }

}
