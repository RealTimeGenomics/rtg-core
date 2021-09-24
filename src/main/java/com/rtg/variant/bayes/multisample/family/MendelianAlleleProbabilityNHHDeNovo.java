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
public final class MendelianAlleleProbabilityNHHDeNovo extends MendelianAlleleProbabilityDeNovo {

  // todo Using log(1/3) for the de novo is only correct for SNPs. It should really do something like log(1/(codesize-1))

  static final double LOG_THIRD = Math.log(1.0 / 3);

  /** Father really is the father. */
  public static final MendelianAlleleProbability SINGLETON_HN = new MendelianAlleleProbabilityNHHDeNovo();

  /** Father really is the mother. */
  public static final MendelianAlleleProbability SINGLETON_NH = new MendelianAlleleProbabilityDeNovo() {
    @Override
    public double probabilityLn(Code code, int father, int mother, int child) {
      return SINGLETON_HN.probabilityLn(code, mother, father, child);
    }
  };

  @Override
  public double probabilityLn(Code code, int father, int mother, int child) {
    assert code.homozygous(father);
    assert code.homozygous(child);
//    assert mother == -1; // Not the case in FB
    return father == child ? Double.NEGATIVE_INFINITY : LOG_THIRD;
  }

}
