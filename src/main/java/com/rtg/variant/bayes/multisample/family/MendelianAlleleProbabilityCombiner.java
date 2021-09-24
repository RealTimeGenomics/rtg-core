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
public final class MendelianAlleleProbabilityCombiner implements MendelianAlleleProbability {

  private final MendelianAlleleProbability mMendelian;
  private final MendelianAlleleProbability mDeNovo;
  private final double mLogRefDenovoPrior;
  private final double mLogNonRefDenovoPrior;
  private final int mReference;

  /**
   * Create a new combined probability M(H) + &#956; * M'(H)
   * @param mendelian singleton for mendelian probability
   * @param denovo singleton for de novo
   * @param logRefDenovoPrior log of prior for de novo
   * @param logNonRefDenovoPrior log of prior for de novo
   * @param ref reference hypothesis
   */
  public MendelianAlleleProbabilityCombiner(MendelianAlleleProbability mendelian, MendelianAlleleProbability denovo, double logRefDenovoPrior, double logNonRefDenovoPrior, int ref) {
    mMendelian = mendelian;
    mDeNovo = denovo;
    mLogRefDenovoPrior = logRefDenovoPrior;
    mLogNonRefDenovoPrior = logNonRefDenovoPrior;
    mReference = ref;
  }

  @Override
  public double probabilityLn(Code code, int father, int mother, int child) {

    // from the docs you either have a valid M(H) or M'(H)
    final double mendelianProbabilityLn = mMendelian.probabilityLn(code, father, mother, child);
    if (mendelianProbabilityLn == Double.NEGATIVE_INFINITY && mLogRefDenovoPrior > Double.NEGATIVE_INFINITY) {
      final double denovo = mDeNovo.probabilityLn(code, father, mother, child);
      final double prior = (father == mother && mother == mReference) ? mLogRefDenovoPrior : mLogNonRefDenovoPrior;
      return prior + denovo;
    }
    return mendelianProbabilityLn;
  }

  /**
   *
   * @param code description of encoding used for nucleotides.
   * @param father the hypothesis from the father's sample.
   * @param mother the hypothesis from the mother's sample. Can be -1 when there is no relevant code (eg mothers contribution to Y chromosome).
   * @param child the hypothesis from a child's sample.
   * @return if the child has de novo allele
   */
  @Override
  public boolean isDenovo(Code code, int father, int mother, int child) {
    if (mLogRefDenovoPrior > Double.NEGATIVE_INFINITY) {
      return mDeNovo.probabilityLn(code, father, mother, child) != Double.NEGATIVE_INFINITY;
    }
    return false;
  }

}
