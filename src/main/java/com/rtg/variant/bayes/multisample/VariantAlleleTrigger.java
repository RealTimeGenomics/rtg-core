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

package com.rtg.variant.bayes.multisample;

import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.Description;

/**
 * Wrapper class for Variant Allele threshold checks
 */
public class VariantAlleleTrigger {
  final double mMinVad;
  final double mMinVaf;

  /**
   * Create a trigger with the specified thresholds.
   * If both thresholds are zero then no check will be made.
   * @param depth Require that variant alleles occur more than {@code depth} (error-adjusted) times. Should be &gt;= 0
   * @param fraction  require that the variant makes up at least this fraction of total evidence. Should be between 0.0 and 1.0.
   */
  public VariantAlleleTrigger(double depth, double fraction) {
    assert fraction >= 0.0 && fraction <= 1.0;
    assert depth >= 0;
    mMinVad = depth;
    mMinVaf = fraction;
  }

  /**
   * Select the highest frequency variant allele if it passes the thresholds
   * @param ac counts and errors of the alleles at this position
   * @param description current description
   * @param refAllele the allele present in the reference sequence
   * @return -1 if there is no variant allele meeting the thresholds or thresholds are disabled. Otherwise the index of the variant allele
   */
  public int getVariantAllele(AlleleStatistics<?> ac, Description description, String refAllele) {
    if (mMinVad <= 0 && mMinVaf <= 0) {
      return -1;
    }
    final int refCode = description.indexOf(refAllele);
    double tot = 0;
    double vac = 0;
    int va = -1;
//    boolean tied = false;
    for (int i = 0; i < description.size(); ++i) {
      final double count = ac.count(i) - ac.error(i);
      tot += count;
      //System.err.println("a = " + i + " (" + description.name(i) + ") ac = " + count);
      if (i != refCode && count >= vac) {
//        tied = count == vac;
        vac = count;
        va = i;
      }
    }
    if (va != -1) {
      final double vaf = vac / tot;
      if (vac < mMinVad || vaf < mMinVaf) {
        va = -1;
      }
    }
//    if (va != -1 && tied) {
//      System.err.println("Two non-ref alleles met VAC/VAF threshold at " + templateName + ":" + (position + 1));
//      System.err.println("va = " + va + " (" + description.name(va) + ") vac = " + vac + " vaf = " + vaf);
//    }
    return va;
  }
}
