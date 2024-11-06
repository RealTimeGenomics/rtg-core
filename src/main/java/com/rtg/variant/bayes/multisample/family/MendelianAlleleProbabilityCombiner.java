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
