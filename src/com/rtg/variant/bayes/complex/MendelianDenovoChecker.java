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
package com.rtg.variant.bayes.complex;

import com.rtg.reference.Ploidy;
import com.rtg.relation.ChildFamilyLookup;
import com.rtg.relation.Family;
import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantSample;

/**
 * Checks for de novo calls in mendelian inheritance scenarios
 */
public class MendelianDenovoChecker implements DenovoChecker {
  final ChildFamilyLookup mLookup;

  /**
   * Constructor
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
    final String[] childAlleles = StringUtils.split(child.getName(), ':');
    //currently VariantSample == null for Ploidy.NONE
    final String[] parentAAlleles = parentA == null || parentA.getPloidy() == Ploidy.NONE ? new String[0] : StringUtils.split(parentA.getName(), ':');
    final String[] parentBAlleles = parentB == null || parentB.getPloidy() == Ploidy.NONE ? new String[0] : StringUtils.split(parentB.getName(), ':');
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
      return     (hasAllele(childCall[0], parentACall) && hasAllele(childCall[1], parentBCall))
                 || (hasAllele(childCall[0], parentBCall) && hasAllele(childCall[1], parentACall));
    }
  }

  private static boolean hasAllele(String s, String[] parentCall) {
    boolean hasAllele = false;
    for (String allele : parentCall) {
      if (s.equals(allele)) {
        hasAllele = true;
      }
    }
    return hasAllele;
  }
}
