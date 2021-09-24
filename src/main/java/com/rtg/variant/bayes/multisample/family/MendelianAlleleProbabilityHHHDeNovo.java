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
 * De novo probabilities for polyploid.
 * father, mother and child treated as haploid
 */
public class MendelianAlleleProbabilityHHHDeNovo extends MendelianAlleleProbabilityDeNovo {

  /** Assuming child inherits from mother */
  public static final MendelianAlleleProbability SINGLETON_HHH = new MendelianAlleleProbabilityHHHDeNovo();

  @Override
  public double probabilityLn(Code code, int father, int mother, int child) {
    assert code.homozygous(father);
    assert code.homozygous(mother);
    assert code.homozygous(child);
    return mother == child ? Double.NEGATIVE_INFINITY : MendelianAlleleProbabilityNHHDeNovo.LOG_THIRD;
  }
}
