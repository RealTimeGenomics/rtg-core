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
import com.rtg.util.integrity.Exam;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.Statistics;

/**
 * Accumulates posteriors for the cancer side of a contaminated cancer sample.
 */
public class ModelCancerContamination extends Model<Description> {

  private final Hypotheses<?> mSubHypotheses;

  private final double mContamination;

  private final double mContaminationM;

  /**
   * @param hyp used to get the probabilities in the various categories for both normal and cancer.
   * @param contamination fraction of cancer sample that is contamination from normal sample.
   * @param statistics for statistics collection
   */
  public ModelCancerContamination(final HypothesesCancer hyp, final double contamination, Statistics<?> statistics) {
    super(hyp, statistics);
    mSubHypotheses = hyp.subHypotheses();
    mContamination = contamination;
    mContaminationM = 1.0 - mContamination;
  }

  private double[] oneDprob(EvidenceInterface evidence) {
    final Code code = mSubHypotheses.code();
    final double[] probs = new double[mSubHypotheses.size()];
    for (int i = 0; i < probs.length; i++) {
      probs[i] = 0.5 * (evidence.probability(code.a(i)) + evidence.probability(code.bc(i)));
    }
    //System.err.println(probs.length + " :" + Arrays.toString(probs));
    return probs;
  }

  @Override
  public void increment(EvidenceInterface evidence) {
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
      final double prob = probs[a] * mContamination + probs[b] * mContaminationM;
      // Phred scores of 0 can result in 0 probabilty, just skip them
      if (prob <= 0.0) {
        return;
      }
      //adjust for mapQ - see theory in scoring.tex
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

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertEquals(1.0, mContamination + mContaminationM, 0.000001);
    Exam.assertTrue(0.0 <= mContamination && mContamination <= 1.0);
    return true;
  }
}
