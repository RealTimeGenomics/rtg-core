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
    final int hypPad = nameLength() + 1;
    final StringBuilder sb = new StringBuilder();
    sb.append("Hypotheses").append(StringUtils.LS);
    for (int i = 0; i < size(); i++) {
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
