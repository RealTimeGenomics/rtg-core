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
 * Mother ignored, father haploid, daughter ignored.
 */
public final class MendelianAlleleProbabilityNHNDeNovo extends MendelianAlleleProbabilityDeNovo {

  /** Father really is the father. */
  public static final MendelianAlleleProbability SINGLETON_HN = new MendelianAlleleProbabilityNHNDeNovo();

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
//    assert mother == -1; //Annonyingly not the case when using FB
//    assert child == -1;
    return Double.NEGATIVE_INFINITY;
  }

}
