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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.util.StringUtils;
import com.rtg.util.format.FormatReal;
import com.rtg.variant.bayes.AlleleBalanceProbability;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.Statistics;

/**
 * Accumulates posteriors for the cancer side of a contaminated cancer sample.
 * @param <S> type of underlying hypotheses
 */
public class ModelCancerContamination<S extends Hypotheses<? extends Description>> extends Model<Description> {
  private final S mSubHypotheses;

  private final double mContamination;

  private final double mContaminationM;

  /**
   * @param hyp used to get the probabilities in the various categories for both normal and cancer.
   * @param contamination fraction of cancer sample that is contamination from normal sample.
   * @param statistics for statistics collection
   * @param alleleBalance allele balance probability calculator
   */
  public ModelCancerContamination(final HypothesesCancer<S> hyp, final double contamination, final Statistics<?> statistics, AlleleBalanceProbability alleleBalance) {
    super(hyp, statistics, alleleBalance);
    if (contamination < 0 || contamination > 1) {
      throw new IllegalArgumentException("Illegal contamination: " + contamination);
    }
    mSubHypotheses = hyp.subHypotheses();
    mContamination = contamination;
    mContaminationM = 1.0 - mContamination;
  }

  private ModelCancerContamination(ModelCancerContamination<S> original) {
    super(original);
    mSubHypotheses = original.mSubHypotheses;
    mContamination = original.mContamination;
    mContaminationM = original.mContaminationM;

  }

  @Override
  public void freeze() {
    // Avoid doing the allele balance calculation of the standard model since it doesn't play
    // well with the cross-product hypotheses used in this model.
    assert !mFrozen : "Should only freeze once";
    mFrozen = true;
  }

  private double[] oneDprob(final EvidenceInterface evidence) {
    final Code code = mSubHypotheses.code();
    final double[] probs = new double[mSubHypotheses.size()];
    for (int i = 0; i < probs.length; ++i) {
      final double a = Math.max(0, evidence.probability(code.a(i)));
      final double b = Math.max(0, evidence.probability(code.bc(i)));
      probs[i] = 0.5 * (a + b);
    }
    return probs;
  }

  @Override
  public void increment(final EvidenceInterface evidence) {
    incrementStatistics(evidence);
    if (ambiguityShortCircuit(evidence)) {
      return;
    }
    final double r = evidence.mapError();
    final double rc = 1.0 - r;
    // Avoid case where mapq is 0 which gives a NaN
    if (rc <= 0.0) {
      return;
    }
    final double[] probs = oneDprob(evidence);
    final double pE = evidence.pe();
    final double pEr = r * pE;
    final Code code = hypotheses().code();
    for (int i = 0; i < size(); ++i) {
      final int a = code.a(i);
      final int b = code.bc(i);
      // The cross-product is normal x cancer
      final double prob = probs[a] * mContamination + probs[b] * mContaminationM;
      // Phred scores of 0 can result in 0 probability, just skip them
      if (prob <= 0.0) {
        return;
      }
      // Adjust for mapQ - see theory in scoring.tex
      final double pr = prob * rc + pEr;
      final double np = arithmetic().multiply(mPosteriors[i], arithmetic().prob2Poss(pr));
      assert arithmetic().isValidPoss(np); //: System.err.println("np=" + np + " mPosteriors[i]=" + mPosteriors[i] + " i=" + i + " pr=" + pr + " prob=" + prob + " rc=" + rc + " r=" + r + " pE=" + pE);
      mPosteriors[i] = np;
    }
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Contaminated Cancer Model");
    final FormatReal fmt = new FormatReal(4, 3);
    sb.append(" contamination=").append(fmt.format(mContamination));
    sb.append(LS);
    final int pad = hypotheses().maxNameLength();
    final int size = ((HypothesesCancer<?>) hypotheses()).subHypotheses().size();
    for (int i = 0; i < size; ++i) {
      sb.append(StringUtils.padLeft(hypotheses().name(i), pad));
      for (int j = 0; j < size; ++j) {
        final int k = hypotheses().code().code(i, j);
        sb.append(fmt.format(arithmetic().poss2Ln(mPosteriors[k])));
      }
      sb.append(LS);
    }
  }

  @Override
  public ModelCancerContamination<S> copy() {
    return new ModelCancerContamination<>(this);
  }
}
