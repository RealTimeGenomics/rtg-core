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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.util.StringUtils;
import com.rtg.util.format.FormatReal;
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
   */
  public ModelCancerContamination(final HypothesesCancer<S> hyp, final double contamination, final Statistics<?> statistics) {
    super(hyp, statistics);
    if (contamination < 0 || contamination > 1) {
      throw new IllegalArgumentException("Illegal contamination: " + contamination);
    }
    mSubHypotheses = hyp.subHypotheses();
    mContamination = contamination;
    mContaminationM = 1.0 - mContamination;
  }

  private double[] oneDprob(final EvidenceInterface evidence) {
    final Code code = mSubHypotheses.code();
    final double[] probs = new double[mSubHypotheses.size()];
    for (int i = 0; i < probs.length; i++) {
      final double a = Math.max(0, evidence.probability(code.a(i)));
      final double b = Math.max(0, evidence.probability(code.bc(i)));
      probs[i] = 0.5 * (a + b);
    }
    //System.err.println(probs.length + " :" + Arrays.toString(probs));
    return probs;
  }

  @Override
  public void increment(final EvidenceInterface evidence) {
    //System.err.println(evidence);
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
    for (int i = 0; i < size(); i++) {
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
    final int pad = hypotheses().nameLength();
    final int size = ((HypothesesCancer) hypotheses()).subHypotheses().size();
    for (int i = 0; i < size; i++) {
      sb.append(StringUtils.padLeft(hypotheses().name(i), pad));
      for (int j = 0; j < size; j++) {
        final int k = hypotheses().code().code(i, j);
        sb.append(fmt.format(arithmetic().poss2Ln(mPosteriors[k])));
      }
      sb.append(LS);
    }
  }

}
