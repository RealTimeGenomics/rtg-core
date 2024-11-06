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
