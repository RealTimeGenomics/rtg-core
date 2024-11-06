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

package com.rtg.variant.bayes.multisample.population;


import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.MultisampleJointScorer;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 */
public class Convergence {

  private final double mErrorRate;

  static int error(double err, int read0, Hypotheses<?> haploid, Random random) {
    final double rand = random.nextDouble();
    if (rand > err) {
      return read0;
    }
    final int n = haploid.size() - 1;
    final int other = random.nextInt(n);
    //System.err.println("read0=" + read0 + " n=" + n + " other=" + other);
    if (other < read0) {
      return other;
    }
    return other + 1;
  }

  static int choose(int i, Hypotheses<?> hyp, Random random) {
    final int a = hyp.code().a(i);
    final int b = hyp.code().bc(i);
    if (a == b) {
      return a;
    }
    final double rand = random.nextDouble();
    if (rand <= 0.5) {
      return a;
    } else {
      return b;
    }
  }

  static int sample(double[] haploid, final Code code, Random random) {
    final int a = sampleHaploid(haploid, random);
    final int b = sampleHaploid(haploid, random);
    return code.code(a, b);
  }

  static int sampleHaploid(double[] distr, Random random) {
    final double rand = random.nextDouble();
    double soFar = 0.0;
    for (int i = 0; i < distr.length; ++i) {
      soFar += distr[i];
      if (soFar >= rand) {
        return i;
      }
    }
    return distr.length - 1;
  }

  static double[] distr(HypothesesPrior<?> hyp) {
    final double[] distr = new double[hyp.size()];
    final PossibilityArithmetic arith = hyp.arithmetic();
    for (int i = 0; i < hyp.size(); ++i) {
      distr[i] = arith.poss2Prob(hyp.p(i));
    }
    return distr;
  }

  /**
   */
  public static class SimulationResult {
    private final int[] mIncorrect;
    private final int mTotalIncorrect;
    /**
     * @param incorrect counts of incorrect calls indexed on hypotheses.
     */
    SimulationResult(int[] incorrect) {
      mIncorrect = incorrect;
      int total = 0;
      for (int i = 0; i < incorrect.length; ++i) {
        total += mIncorrect[i];
      }
      mTotalIncorrect = total;
    }
    /**
     * Get totalIncorrect.
     * @return Returns the totalIncorrect.
     */
    int totalIncorrect() {
      return mTotalIncorrect;
    }

    int[] incorrect() {
      return mIncorrect;
    }
  }

  private final HypothesesPrior<Description> mDiploid;

  private final HypothesesPrior<Description> mHaploid;

  private final Description mDescription;

  private final List<ModelInterface<?>> mModels;

  private final Random mRandom;

  private final int[] mSample;

  private final MultisampleJointScorer mPopulation;

  /**
   * Case when actual distribution is equal to the priors.
   * @param haploid set of potential hypotheses.
   * @param diploid set of potential hypotheses.
   * @param estimator given the actual calls re-estimates the parameters.
   * @param models which are updated and make calls.
   * @param random number generator.
   */
  Convergence(final HypothesesPrior<Description> haploid, final HypothesesPrior<Description> diploid, final Estimator estimator, final List<ModelInterface<?>> models, final Random random) {
    this(haploid, diploid, distr(haploid), 0.05, estimator, models, random);
  }

  /**
   * @param haploid set of potential hypotheses.
   * @param diploid set of potential hypotheses.
   * @param distr actual distribution of (haploid) probabilities (may differ from th priors in the hypotheses).
   * @param errorRate error rate for individual nucleotides.
   * @param estimator given the actual calls re-estimates the parameters.
   * @param models which are updated and make calls.
   * @param random number generator.
   */
  Convergence(final HypothesesPrior<Description> haploid, final HypothesesPrior<Description> diploid, final double[] distr, final double errorRate, final Estimator estimator, final List<ModelInterface<?>> models, final Random random) {
    mDiploid = diploid;
    mHaploid = haploid;
    mDescription = mHaploid.description();
    mErrorRate = errorRate;
    mModels = models;
    mRandom = random;
    mSample = new int[models.size()];
    for (int i = 0; i < models.size(); ++i) {
      final ModelInterface<?> model = models.get(i);
      final int sample;
      if (model.haploid()) {
        sample = sampleHaploid(distr, mRandom);
      } else {
        sample = sample(distr, mDiploid.code(), mRandom);
      }
      mSample[i] = sample;
    }
    mPopulation = new EmAlgorithm(estimator, 50);
  }

  SimulationResult simulate() {
    for (int i = 0; i < mModels.size(); ++i) {
      final ModelInterface<?> model = mModels.get(i);
      final Hypotheses<?> hyp = model.hypotheses();
      final int read0 = choose(mSample[i], hyp, mRandom);
      final int read = error(0.01, read0, mHaploid, mRandom);
      final EvidenceInterface ev = new EvidenceQ(mDescription, read, 0, 0, 0.01, mErrorRate, true, false, true, false, false);
      model.increment(ev);
    }
    final List<ModelInterface<?>> models = mModels.stream().map(ModelInterface::copy).collect(Collectors.toList());
    models.stream().forEach(ModelInterface::freeze);
    final HypothesisScores popCalls = mPopulation.getBestScores(models, new PriorContainer<>(new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, mHaploid, mDiploid), null));
    final int[] incorrect = new int[mDiploid.size()];
    for (int i = 0; i < mModels.size(); ++i) {
      final int sample = mSample[i];
      if (popCalls.getScores()[i].hypothesis() != sample) {
        incorrect[sample]++;
      }
    }
    return new SimulationResult(incorrect);
  }

}
