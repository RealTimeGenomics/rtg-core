/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.HypothesesPowerSet;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * An allele based somatic caller.
 */
public class AlleleSomaticCaller extends AbstractSomaticCaller {

  /**
   * @param qFactory Q matrix factory
   * @param params variant params
   * @param phi probability of seeing contrary evidence in the original
   * @param psi probability of seeing contrary evidence in the derived
   */
  public AlleleSomaticCaller(AlleleSomaticPriorsFactory<?> qFactory, VariantParams params, double phi, double psi) {
    super(qFactory, qFactory, params, phi, psi);
  }

  @Override
  protected AbstractSomaticPosterior makePosterior(final ModelInterface<?> normal, final ModelInterface<?> cancer, HypothesesPrior<?> hypotheses, double mu) {
    return new SomaticPosteriorAllele(normal.hypotheses(), (HypothesesPowerSet<?>) cancer.hypotheses(), mPhi, mPsi);
  }
}
