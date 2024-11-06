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
package com.rtg.variant.bayes.complex;

import com.rtg.reference.Ploidy;
import com.rtg.relation.ChildFamilyLookup;
import com.rtg.relation.Family;
import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantSample;
import com.rtg.variant.util.VariantUtils;

/**
 * Checks for de novo calls in mendelian inheritance scenarios.
 */
public class MendelianDenovoChecker implements DenovoChecker {
  final ChildFamilyLookup mLookup;

  /**
   * Constructor
   *
   * @param lookup family lookup table
   */
  public MendelianDenovoChecker(ChildFamilyLookup lookup) {
    mLookup = lookup;
  }


  @Override
  public boolean isDenovo(Variant variant, int sample) {
    final VariantSample child = variant.getSample(sample);
    final Family family = mLookup.getFamily(sample);
    final VariantSample parentA = variant.getSample(family.getSampleIds()[Family.FATHER_INDEX]);
    final VariantSample parentB = variant.getSample(family.getSampleIds()[Family.MOTHER_INDEX]);
    final String[] childAlleles = StringUtils.split(child.getName(), VariantUtils.COLON);
    //currently VariantSample == null for Ploidy.NONE
    final String[] parentAAlleles = parentA == null || parentA.getPloidy() == Ploidy.NONE ? new String[0] : StringUtils.split(parentA.getName(), VariantUtils.COLON);
    final String[] parentBAlleles = parentB == null || parentB.getPloidy() == Ploidy.NONE ? new String[0] : StringUtils.split(parentB.getName(), VariantUtils.COLON);
    if (child.getPloidy() == Ploidy.DIPLOID) {
      return !isMendelianDiploid(childAlleles, parentAAlleles, parentBAlleles);
    } else if (child.getPloidy() == Ploidy.HAPLOID) {
      return !isMendelianHaploid(childAlleles, parentAAlleles, parentBAlleles);
    } else {
      throw new UnsupportedOperationException("Can only handle haploid and diploid");
    }
  }

  private static boolean isMendelianHaploid(String[] childAlleles, String[] parentAAlleles, String[] parentBAlleles) {
    if (parentAAlleles.length > 1) {
      return hasAllele(childAlleles[0], parentAAlleles);
    } else if (parentBAlleles.length > 1) {
      return hasAllele(childAlleles[0], parentBAlleles);
    } else {
      return hasAllele(childAlleles[0], parentAAlleles) || hasAllele(childAlleles[0], parentBAlleles);
    }
  }

  private static boolean isMendelianDiploid(String[] childCall, String[] parentACall, String[] parentBCall) {
    assert childCall.length == 2;
    if (childCall[0].equals(childCall[1])) {
      //homozygous call, both parents must have allele
      return hasAllele(childCall[0], parentACall) && hasAllele(childCall[0], parentBCall);
    } else {
      //one allele from each parent
      return (hasAllele(childCall[0], parentACall) && hasAllele(childCall[1], parentBCall))
        || (hasAllele(childCall[0], parentBCall) && hasAllele(childCall[1], parentACall));
    }
  }

  private static boolean hasAllele(String s, String[] parentCall) {
    boolean hasAllele = false;
    for (final String allele : parentCall) {
      if (s.equals(allele)) {
        hasAllele = true;
      }
    }
    return hasAllele;
  }
}
