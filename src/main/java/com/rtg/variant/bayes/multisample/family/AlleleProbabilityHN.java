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
 */
public final class AlleleProbabilityHN implements AlleleProbability {

  private static final double LN_HALF = Math.log(0.5);

  /** First parent is haploid, second is empty. */
  public static final AlleleProbability SINGLETON_HN = new AlleleProbabilityHN();

  /** Second parent is haploid, first is empty. */
  public static final AlleleProbability SINGLETON_NH = new AlleleProbability() {
    @Override
    public double probabilityLn(Code code, int ref, int father, int mother) {
      return SINGLETON_HN.probabilityLn(code, ref, mother, father);
    }
  };

  @Override
  public double probabilityLn(Code code, int ref, int father, int mother) {
    assert code.homozygous(father);
//    assert mother == -1; // Not the case in FB
    return ref == father || ref == -1 ? 0.0 : LN_HALF; //log(1.0)
  }

}
