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


/**
 * @param <D> description type
 * @param <H> class of <code>Hypotheses</code>. Needed so <code>ModelCancerFactory</code> can insist on a <code>HypothesesCancer</code>.
 */
public interface ModelFactory<D extends Description, H extends Hypotheses<D>> extends ReferenceBasedFactory<ModelInterface<D>> {

  /**
   * Make a model given the code for the reference sequence.
   *
   * @param ref the reference code.
   * @return the model.
   */
  @Override
  ModelInterface<D> make(int ref);

  /**
   * Return the same hypotheses as used in a model returned by make.
   * Note that the priors contained in these hypotheses should not generally be used
   * (as they do not take into account site-specific priors, EM iterations etc).
   * @param ref the reference code.
   * @return the hypotheses.
   */
  H defaultHypotheses(int ref);
}
