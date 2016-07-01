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

  private final AlleleBalanceProbability mAlleleBalance;

  protected boolean mFrozen = false;

  protected boolean ambiguityShortCircuit(final EvidenceInterface evidence) {
    return evidence.mapError() >= AMBIGUITY_THRESHOLD; // || evidence.error() > MAX_BASE_ERROR;
  }

  /**
   * @param hyp provides information about the set of hypotheses.
   * @param statistics to store counts
   * @param alleleBalance allele balance probability implementation
   */
  public Model(final Hypotheses<D> hyp, final Statistics<?> statistics, AlleleBalanceProbability alleleBalance) {
    mHypotheses = hyp;
    mStatistics = statistics;
    final int size = hyp.code().size();
    mPosteriors = new double[size];
    Arrays.fill(mPosteriors, arithmetic().one());
    mAlleleBalance = alleleBalance;
    //    System.err.println("construct" + hypotheses.getClass().getName());
    //    System.err.println("hyparith" + hypotheses.arithmetic().getClass().getName());
  }

  /**
   * test copy constructor
   * @param m model to copy
   */
  protected Model(Model<D> m) {
    mHypotheses = m.mHypotheses;
    mStatistics = (Statistics<?>) m.mStatistics.copy();
    mPosteriors = Arrays.copyOf(m.mPosteriors, m.mPosteriors.length);
    mAlleleBalance = m.mAlleleBalance;
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
    assert !mFrozen : "This model has been frozen you should not be updating it any more";
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

  @Override
  public void freeze() {
    assert !mFrozen : "Should only freeze once";
    mFrozen = true;
    final PossibilityArithmetic arithmetic = arithmetic();
    for (int hyp = 0; hyp < mPosteriors.length; hyp++) {
      mPosteriors[hyp] = arithmetic.multiply(mPosteriors[hyp], arithmetic.ln2Poss(mAlleleBalance.alleleBalanceLn(hyp, hypotheses(), statistics())));
    }
  }

  protected void incrementStatistics(EvidenceInterface distribution) {
    mStatistics.increment(distribution, mHypotheses.reference());
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

  @Override
  public AlleleBalanceProbability alleleBalanceProbability() {
    return mAlleleBalance;
  }

  private double posterior(final PossibilityArithmetic arith, final int i, final Factor<?> hypotheses) {
    assert mFrozen : "You should freeze the model before calling posterior methods";
    return arith.multiply(mPosteriors[i], hypotheses.p(i));
  }

  @Override
  public double posteriorLn0(int hyp) {
    assert mFrozen : "You should freeze the model before calling posterior methods";
    return arithmetic().poss2Ln(mPosteriors[hyp]);
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

  @Override
  public Model<D> copy() {
    return new Model<>(this);
  }

}
