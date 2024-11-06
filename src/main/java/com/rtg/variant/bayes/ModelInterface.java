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
