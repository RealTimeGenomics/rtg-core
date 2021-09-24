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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.variant.bayes.Code;

/**
 * Cases for the Y chromosome in mammals.
 * Mother ignored, father and child (son) haploid.
 */
public final class MendelianAlleleProbabilityNHH extends MendelianAlleleProbabilityNormal {

  /** Father really is the father. */
  public static final MendelianAlleleProbability SINGLETON_HN = new MendelianAlleleProbabilityNHH();

  /** Father really is the mother. */
  public static final MendelianAlleleProbability SINGLETON_NH = new MendelianAlleleProbabilityNormal() {
    @Override
    public double probabilityLn(Code code, int father, int mother, int child) {
      return SINGLETON_HN.probabilityLn(code, mother, father, child);
    }
  };

  @Override
  public double probabilityLn(Code code, int father, int mother, int child) {
    assert code.homozygous(father);
    assert code.homozygous(child);
    //father and mother have been switched around in the NHH case (see SINGLETON_NH)
 //   assert mother == -1; //this is prohibitively annoying to ensure when doing forward backward
    return father == child ? 0.0 : Double.NEGATIVE_INFINITY;
  }

}
