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


import java.util.Arrays;

import com.rtg.launcher.GlobalFlags;
import com.rtg.reference.Ploidy;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * A Bayesian model of a set of hypotheses selected by an integer index.
 * @param <D> description type
 */
public class Model<D extends Description> extends IntegralAbstract implements ModelInterface<D> {

  // Factor 1.1 to cover arithmetic error in sum used for error() in EvidenceComplex
  private static final double MAX_BASE_ERROR = VariantUtils.phredToProb(GlobalFlags.getIntegerValue(GlobalFlags.MIN_BASE_QUALITY)) * 1.1;
  private static final double EXPECTED_ALLELE_FREQUENCY = GlobalFlags.getDoubleValue(GlobalFlags.CALLER_EXPECTED_ALLELE_FREQUENCY);

  /**
   * When <code>MAPQ</code> is less than or equal to this then treat as ambiguous.
   */
  public static final int AMBIGUITY_PHRED = 3;

  /**
   * When <code>mapError()</code> is greater than or equal to this then treat as ambiguous.
   */
  public static final double AMBIGUITY_THRESHOLD = VariantUtils.phredToProb(AMBIGUITY_PHRED);

  private static final String SPACES3 = "   ";

  private static final String SPACES4 = "    ";

  private final Hypotheses<D> mHypotheses;

  protected final double[] mPosteriors;

  private final Statistics<?> mStatistics;

  protected boolean ambiguityShortCircuit(final EvidenceInterface evidence) {
    return evidence.mapError() >= AMBIGUITY_THRESHOLD || evidence.error() > MAX_BASE_ERROR;
  }

  /**
   * @param hyp provides information about the set of hypotheses.
   * @param statistics to store counts
   */
  public Model(final Hypotheses<D> hyp, final Statistics<?> statistics) {
    mHypotheses = hyp;
    mStatistics = statistics;
    final int size = hyp.code().size();
    mPosteriors = new double[size];
    Arrays.fill(mPosteriors, arithmetic().one());
    //    System.err.println("construct" + hypotheses.getClass().getName());
    //    System.err.println("hyparith" + hypotheses.arithmetic().getClass().getName());
  }

  @Override
  public final Hypotheses<D> hypotheses() {
    return mHypotheses;
  }

  @Override
  public final int reference() {
    return mHypotheses.reference();
  }

  @Override
  public final int size() {
    return mHypotheses.code().size();
  }

  @Override
  public final String name(int i) {
    return mHypotheses.name(i);
  }

  @Override
  public final Description description() {
    return mHypotheses.description();
  }

  @Override
  public final boolean haploid() {
    return mHypotheses.haploid();
  }

  /**
   * @return the possibility arithmetic used in this model.
   */
  @Override
  public final PossibilityArithmetic arithmetic() {
    return mHypotheses.arithmetic();
  }

  private void increment(EvidenceQ evidence) {
    for (int i = 0; i < size(); i++) {
      final double v = evidence.logEvidentialProbability(i);
      if (v > 0) {
        return;
      }
      mPosteriors[i] = arithmetic().multiply(mPosteriors[i], arithmetic().ln2Poss(v));
    }
  }

  @Override
  public void increment(EvidenceInterface evidence) {
    incrementStatistics(evidence);
    if (ambiguityShortCircuit(evidence)) {
      return;
    }
    if (evidence instanceof EvidenceQ) {
      // EvidenceQ objects have precomputed evidence probabilities
      increment((EvidenceQ) evidence);
    } else {
      final double r = evidence.mapError();
      final double rc = 1.0 - r;
      // Avoid case where mapq is 0 which gives a NaN
      if (rc <= 0.0) {
        return;
      }
      final Code code = hypotheses().code();
      final double pE = evidence.pe();
      //assert pE >= 0 && pE <= 1;
      final double pEr = r * pE;
      for (int i = 0; i < size(); i++) {
        final double prob = 0.5 * (evidence.probability(code.a(i)) + evidence.probability(code.bc(i)));
        // Phred scores of 0 can result in 0 probabilty, just skip them
        if (prob <= 0.0) {
          return;
        }

        final double pr = prob * rc + pEr;
        //System.err.println("i=" + i + " pr=" + pr + " prob=" + prob + " rc=" + rc + " r=" + r + " pE=" + pE);
        //assert pr > 0.0; // : "pr=" + pr + " prob=" + prob + " rc=" + rc + " r=" + r + " pE=" + pE;
        final double np = arithmetic().multiply(mPosteriors[i], arithmetic().prob2Poss(pr));
        //assert arithmetic().isValidPoss(np);
        mPosteriors[i] = np;
      }
    }
  }

  protected void incrementStatistics(EvidenceInterface distribution) {
    mStatistics.increment(distribution, mHypotheses.reference());
  }

  protected double alleleBalanceLn(int i) {
    if (/* hypotheses().ploidy() != Ploidy.DIPLOID && */ hypotheses().ploidy() != Ploidy.HAPLOID) {
      return 0.0;
    }
    final int trials = statistics().coverage();
    if (trials == 0) {
      return 1.0;
    }
    final double abp;
    final int a = hypotheses().code().a(i);
    final int b = hypotheses().code().bc(i);
    assert a == b; // haploid
    final AlleleStatistics<?> counts = statistics().counts();
//    if (a == b) {
    final double p = a == hypotheses().reference() ? 1 - EXPECTED_ALLELE_FREQUENCY : EXPECTED_ALLELE_FREQUENCY;
    //abp = MathUtils.hoeffdingLn(trials, MathUtils.round(counts.count(a)), p);
    abp = -MathUtils.logBinomial(p, trials, (int) MathUtils.round(counts.count(a)));
//    } else {
//      final long observed0 = MathUtils.round(counts.count(a));
//      final long observed1 = MathUtils.round(counts.count(b));
//      abp = MathUtils.hoeffdingLn(trials, observed0, 0.5) + MathUtils.hoeffdingLn(trials, observed1, 0.5);
//    }
    return abp;
  }

  @Override
  public HypothesisScore best(final HypothesesPrior<?> hypotheses) {
    if (hypotheses.size() != size()) {
      throw new RuntimeException("" + " hypotheses.size()=" + hypotheses.size() + " size()=" + size());
    }
    final PossibilityArithmetic arith = arithmetic();
    final double[] posteriors = new double[size()];
    for (int i = 0; i < size(); i++) {
      posteriors[i] = posterior(arith, i, hypotheses);
    }
    return new HypothesisScore(new ArrayGenotypeMeasure(arith, posteriors, hypotheses));
  }

  private double posterior(final PossibilityArithmetic arith, final int i, final Factor<?> hypotheses) {
    return arith.multiply(arith.multiply(mPosteriors[i], arith.ln2Poss(alleleBalanceLn(i))), hypotheses.p(i));
  }

  @Override
  public double posteriorLn0(int hyp) {
    return arithmetic().multiply(mPosteriors[hyp], arithmetic().ln2Poss(alleleBalanceLn(hyp)));
  }

  @Override
  public Statistics<?> statistics() {
    return mStatistics;
  }

  @Override
  public void statistics(StringBuilder sb, final HypothesesPrior<?> hypotheses) {
    final int size = hypotheses().code().size();
    for (int i = 0; i < size; i++) {
      statistics(sb, i, hypotheses);
      sb.append(StringUtils.LS);
    }
  }

  void statistics(StringBuilder sb, final int hyp, final HypothesesPrior<?> hypotheses) {
    assert mHypotheses.code().valid(hyp);
    sb.append(SPACES4);
    sb.append(mHypotheses.reference() == hyp ? '*' : ' ');
    mHypotheses.writeName(sb, hyp);
    sb.append(SPACES3);
    sb.append(Utils.realFormat(posteriorLn0(hyp), 3));
    sb.append(SPACES3);
    sb.append(Utils.realFormat(posterior(arithmetic(), hyp, hypotheses), 3));
  }

  @Override
  public double p(int code) {
    return mPosteriors[code];
  }

  @Override
  public boolean isNormalized() {
    return true;
  }

  @Override
  public Factor<D> normalize() {
    return this;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("Model");
    //mStatistics.toString(sb);
    sb.append(LS);
    final int hypPad = mHypotheses.nameLength();
    for (int i = 0; i < size(); i++) {
      sb.append(StringUtils.padLeft(mHypotheses.name(i), hypPad));
      sb.append(" ").append(StringUtils.padLeft(Utils.realFormat(posteriorLn0(i), 3), 7));
      sb.append(LS);
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mHypotheses);
    return true;
  }

}
