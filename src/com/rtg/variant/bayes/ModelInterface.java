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

package com.rtg.variant.bayes;

import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * A Bayesian model of a set of hypotheses selected by an integer index.
 * @param <D> description type
 */
public interface ModelInterface<D extends Description> extends EvidenceAcceptor, Factor<D> {

  /**
   * Get whether the model is haploid
   * @return true if the model is haploid
   */
  boolean haploid();

  /**
   * Get the hypothesis id corresponding to the reference.
   * @return the reference hypotheses
   */
  int reference();

  /**
   * Get the name of a particular hypothesis
   * @param i the hypothesis
   * @return the name of the hypotheses
   */
  String name(int i);

  /**
   * Get the Description for the hypotheses
   * @return the Description
   */
  Description description();

  /**
   * Create a detailed human readable description of each hypothesis and its posterior.
   *
   * @param sb string builder to write statistics to.
   * @param hypotheses with priors for the calling.
   */
  void statistics(StringBuilder sb, HypothesesPrior<?> hypotheses);

  /**
   * @return statistics accumulated across the increment calls.
   */
  Statistics<?> statistics();

  /**
   * @param hyp the hypothesis.
   * @return the natural log of the posterior (not including prior).
   */
  double posteriorLn0(int hyp);

  /**
   * Call using priors in supplied hypotheses.
   * @param hypotheses with priors for the calling.
   * @return the best scoring call.
   */
  HypothesisScore best(HypothesesPrior<?> hypotheses);

  /**
   * Return the object for computing the allele balance probability used with this model.
   * @return allele balance computer
   */
  AlleleBalanceProbability alleleBalanceProbability();

  /**
   * @return a copy of this model
   */
  ModelInterface<D> copy();
}
