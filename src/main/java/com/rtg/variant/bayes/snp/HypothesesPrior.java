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

package com.rtg.variant.bayes.snp;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Factor;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Provides hypotheses with fixed, array-backed prior probabilities (stored as possibilities).
 * @param <D> class of description.
 */
@TestClass("com.rtg.variant.bayes.snp.HypothesesSnpTest")
public class HypothesesPrior<D extends Description> extends Hypotheses<D> implements Factor<D> {

  protected final int mReference;
  protected final double[] mPriors;

  /**
   * @param description of the haploid
   * @param arithmetic used for reporting priors.
   * @param haploid true iff to be haploid.
   * @param ref the reference hypothesis.
   */
  protected HypothesesPrior(D description, PossibilityArithmetic arithmetic, boolean haploid, int ref) {
    super(description, arithmetic, haploid);
    mPriors = new double[size()];
    mReference = ref;
  }

  /**
   * @param description of the haploid
   * @param arithmetic used for reporting priors.
   * @param priors in possibility format.
   * @param haploid true iff to be haploid.
   * @param ref the reference hypothesis.
   */
  public HypothesesPrior(D description, PossibilityArithmetic arithmetic, double[] priors, boolean haploid, int ref) {
    super(description, arithmetic, haploid);
    mPriors = priors;
    mReference = ref;
  }

  /**
   * @param description of the haploid
   * @param arithmetic used for reporting priors.
   * @param priors in possibility format.
   * @param haploid true iff to be haploid.
   * @param ref the reference hypothesis.
   */
  public HypothesesPrior(D description, PossibilityArithmetic arithmetic, HypothesesPrior<D> priors, boolean haploid, int ref) {
    super(description, arithmetic, haploid);
    mPriors = priors.mPriors;
    mReference = ref;
  }

  @Override
  public int reference() {
    return mReference;
  }

  /**
   * Gets the prior for an hypothesis.
   * @param hyp selects the hypothesis.
   * @return the prior selected by <code>hyp</code>, this will be in the format of <code>arithmetic()</code>.
   */
  @Override
  public double p(int hyp) {
    return mPriors[hyp];
  }

  /**
   * Sets the prior for an hypothesis, should only be called during initialization.
   * @param hyp selects the hypothesis.
   * @param p prior selected for <code>hyp</code>, this will be in the format of <code>arithmetic()</code>.
   */
  protected final void setPrior(int hyp, double p) {
    mPriors[hyp] = p;
  }

  @Override
  public String toString() {
    final int hypPad = maxNameLength() + 1;
    final StringBuilder sb = new StringBuilder();
    sb.append("Hypotheses").append(StringUtils.LS);
    for (int i = 0; i < size(); ++i) {
      sb.append(StringUtils.padLeft(name(i), hypPad));
      final double prior = arithmetic().poss2Ln(p(i));
      sb.append(" ").append(StringUtils.padLeft(Utils.realFormat(prior, 3), 7));
      sb.append(StringUtils.LS);
    }
    return sb.toString();
  }

  @Override
  public Hypotheses<D> hypotheses() {
    return this;
  }

  @Override
  public boolean isNormalized() {
    return true;
  }

  @Override
  public Factor<D> normalize() {
    return this;
  }

}
