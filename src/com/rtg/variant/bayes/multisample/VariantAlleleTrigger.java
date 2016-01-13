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

package com.rtg.variant.bayes.multisample;

import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.Description;

/**
 * Wrapper class for Variant Allele thresholding.
 */
public class VariantAlleleTrigger {
  final int mMinVac;
  final double mMinVaf;

  /**
   * Create a trigger with the specified thresholds.
   * If both thresholds are zero then Variant Allele handling will be skipped
   * @param count Require that variant alleles occur more than {@code count} times. Should be >= 0
   * @param fraction  require that the variant makes up at least this fraction of total evidence. Should be between 0.0 and 0.1.
   */
  public VariantAlleleTrigger(int count, double fraction) {
    assert fraction >= 0.0 && fraction <= 1.0;
    assert count >= 0;
    mMinVac = count;
    mMinVaf = fraction;
  }

  /**
   * Select the highest frequency variant allele if it passes the thresholds
   * @param ac counts and errors of the alleles at this position
   * @param description current description
   * @param refAllele the allele present in the reference sequence
   * @return -1 if there is no variant allele meeting the thresholds. Otherwise the index of the variant allele
   */
  public int getVariantAllele(AlleleStatistics<?> ac, Description description, String refAllele) {
    if (mMinVac <= 0 && mMinVaf <= 0) {
      return -1;
    }
    final int refCode = description.indexOf(refAllele);
    double tot = 0;
    double vac = 0;
    int va = -1;
//    boolean tied = false;
    for (int i = 0; i < description.size(); i++) {
      final double count = ac.count(i) - ac.error(i);
      tot += count;
      //System.err.println("a = " + i + " (" + description.name(i) + ") ac = " + count);
      if (i != refCode && count >= vac) {
//        tied = count == vac;
        vac = count;
        va = i;
      }
    }
    double vaf;
    if (va != -1) {
      vaf = vac / tot;
      if (vac <= mMinVac || vaf <= mMinVaf) {
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