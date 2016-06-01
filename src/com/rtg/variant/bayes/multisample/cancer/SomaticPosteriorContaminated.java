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

import com.rtg.launcher.GlobalFlags;
import com.rtg.util.ChiSquared;
import com.rtg.util.MathUtils;
import com.rtg.variant.bayes.AlleleStatistics;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.Statistics;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Calculate the joint posterior and marginal distributions for cancer.
 */
class SomaticPosteriorContaminated extends AbstractSomaticPosterior {

  /** Select Dirichlet distribution for tumor side. */
  public static final String DIRICHLET = "dirichlet";
  /** Select binomial distribution for tumor side. */
  public static final String BINOMIAL = "binomial";
  /** Select a band to allow. */
  public static final String BAND = "band";
  /** Select a band to allow. */
  public static final String MULTIPLE_BINOMIAL = "multiple-binomial";
  /** You can also pick anything else for no allele balance. */

  private final double mLowerBand;
  private final double mUpperBand;
  private final double mU;
  private final double mV;
  private final double[] mPivots;

  {
    final String abType = GlobalFlags.getStringValue(GlobalFlags.TUMOR_ALLELE_BALANCE);
    if (abType.startsWith(BAND)) {
      final String[] parts = abType.split(",");
      final double a = Double.parseDouble(parts[1]);
      final double b = Double.parseDouble(parts[2]);
      final double w = Double.parseDouble(parts[3]);
      mLowerBand = a;
      mUpperBand = b;
      mU = Math.log(w / ((w - 1) * (b - a) + 1));
      mV = Math.log(1.0 / ((w - 1) * (b - a) + 1));
      mPivots = null;
    } else {
      mLowerBand = 0;
      mUpperBand = 0;
      mU = 0;
      mV = 0;
      if (abType.startsWith(MULTIPLE_BINOMIAL)) {
        final String[] parts = abType.split(",");
        mPivots = new double[parts.length - 1];
        for (int k = 1; k < parts.length; k++) {
          mPivots[k - 1] = Double.parseDouble(parts[k]);
        }
      } else {
        mPivots = null;
      }
    }
  }

  private static double alleleBalanceDirichletProbabilityLn(final Hypotheses<?> hyp, final Statistics<?> statistics, final int normalHyp, final int cancerHyp, final double alpha) {
    final Code code = hyp.code();
    final int na = code.a(normalHyp);
    final int nb = code.bc(normalHyp);
    final int ca = code.a(cancerHyp);
    final int cb = code.bc(cancerHyp);
    final double[] p = new double[hyp.description().size()];
    p[na] += 0.5 * alpha;
    p[nb] += 0.5 * alpha;
    p[ca] += 0.5 * (1 - alpha);
    p[cb] += 0.5 * (1 - alpha);
    double lnP = 0; // i.e. p = 1
    double dp = 0;
    final AlleleStatistics<?> counts = statistics.counts();
    for (int k = 0; k < p.length; k++) {
      if (p[k] > 0) { // Avoid infinity arising from p = 0 situation
        final double c = counts.count(k) - counts.error(k);
        if (c > 0) {
          dp += c;
          lnP += Math.log(p[k]) * (c - 1) - ChiSquared.lgamma(c);
        }
      }
    }
    if (dp > 0) {
      lnP += ChiSquared.lgamma(dp);
    }
    return lnP;
  }

  private static double alleleBalanceMultinomialProbabilityLn(final Hypotheses<?> hyp, final Statistics<?> statistics, final int normalHyp, final int cancerHyp, final double alpha) {
    final Code code = hyp.code();
    final int na = code.a(normalHyp);
    final int nb = code.bc(normalHyp);
    final int ca = code.a(cancerHyp);
    final int cb = code.bc(cancerHyp);
    final double[] p = new double[hyp.description().size()];
    p[na] += 0.5 * alpha;
    p[nb] += 0.5 * alpha;
    p[ca] += 0.5 * (1 - alpha);
    p[cb] += 0.5 * (1 - alpha);
    double lnP = 0; // i.e. p = 1
    double dp = 0;
    final AlleleStatistics<?> counts = statistics.counts();
    for (int k = 0; k < p.length; k++) {
      if (p[k] > 0) { // Avoid infinity arising from p = 0 situation
        final double c = counts.count(k) - counts.error(k);
        dp += c;
        lnP += Math.log(p[k]) * c - ChiSquared.lgamma(c + 1);
      }
    }
    if (dp > 0) {
      lnP += ChiSquared.lgamma(dp + 1);
    }
    return lnP;
  }

  private double alleleBalanceMultinomialPivotsLn(final Hypotheses<?> hyp, final Statistics<?> statistics, final int normalHyp, final int cancerHyp) {
    final Code code = hyp.code();
    final int na = code.a(normalHyp);
    final int nb = code.bc(normalHyp);
    final int ca = code.a(cancerHyp);
    final int cb = code.bc(cancerHyp);
    double lnP = 0; // i.e. p = 1
    final double[] p = new double[hyp.description().size()];
    for (final double alpha : mPivots) {
      p[na] += 0.5 * alpha;
      p[nb] += 0.5 * alpha;
      p[ca] += 0.5 * (1 - alpha);
      p[cb] += 0.5 * (1 - alpha);
      double dp = 0;
      final AlleleStatistics<?> counts = statistics.counts();
      for (int k = 0; k < p.length; k++) {
        if (p[k] > 0) { // Avoid infinity arising from p = 0 situation
          final double c = counts.count(k) - counts.error(k);
          dp += c;
          lnP += Math.log(p[k]) * c - ChiSquared.lgamma(c + 1);
        }
      }
      if (dp > 0) {
        lnP += ChiSquared.lgamma(dp + 1);
      }
    }
    return lnP;
  }

  private double alleleBalanceBandLn(final Hypotheses<?> hyp, final Statistics<?> statistics, final int normalHyp, final int cancerHyp) {
    final Code code = hyp.code();
    final int na = code.a(normalHyp);
    final int nb = code.bc(normalHyp);
    final int ca = code.a(cancerHyp);
    final int cb = code.bc(cancerHyp);
//    final int ref = hyp.reference();
//    final int ra = code.a(ref);
//    final int rb = code.bc(ref);
    final AlleleStatistics<?> counts = statistics.counts();

    final double normalCount = counts.count(na) - counts.error(na) + counts.count(nb) - counts.error(nb);
    double somaticCount = 0;
    boolean isSomaticHyp = false;
    if (ca != na && ca != nb) {
      somaticCount += counts.count(ca) - counts.error(ca);
      isSomaticHyp = true;
    }
    if (cb != na && cb != nb) {
      somaticCount += counts.count(cb) - counts.error(cb);
      isSomaticHyp = true;
    }
    if (!isSomaticHyp) {
      return 0; // Not a somatic hypothesis
    }
    final double somaticFraction = (somaticCount + 0.5) / (normalCount + somaticCount + 1.0); // Laplace
    return somaticFraction >= mLowerBand && somaticFraction < mUpperBand ? mU : mV;
  }

  /**
   * @param qa The Q matrix of cancer mutation probabilities.
   * @param normal model.
   * @param cancer model, has cross product hypotheses to allow for contamination.
   * @param hypotheses the hypotheses containing priors
   * @param phi probability of seeing contrary evidence in the original
   * @param psi probability of seeing contrary evidence in the derived
   * @param alpha contamination
   * @param useAlleleBalanceCorrection true if the expected allele balance correction is to be applied
   */
  SomaticPosteriorContaminated(final double[][] qa, final ModelInterface<?> normal, final ModelInterface<?> cancer, HypothesesPrior<?> hypotheses, double phi, double psi, double alpha, boolean useAlleleBalanceCorrection) {
    super(normal.hypotheses(), phi, psi);
    //System.err.println("normal " + normal);
    //System.err.println("cancer " + cancer);
    assert cancer instanceof ModelCancerContamination;
    assert !cancer.haploid(); // The cancer model is a cross-product of normal x cancer hypotheses
    final Code code = cancer.hypotheses().code();
    for (int i = 0; i < mPosterior.length; i++) {
      final double pi = hypotheses.arithmetic().poss2Ln(hypotheses.p(i)) + normal.posteriorLn0(i);
      for (int j = 0; j < mPosterior.length; j++) {
        final int k = code.code(i, j); // code point for normal(i) x cancer(j)
        assert k >= 0 && k < code.size() : k + " " + code.size();
        final double pj = cancer.posteriorLn0(k);
        final double q = MathUtils.log(qa[i][j]);
        double t = q + pi + pj;
        if (useAlleleBalanceCorrection) {
          final String abType = GlobalFlags.getStringValue(GlobalFlags.TUMOR_ALLELE_BALANCE);
          if (DIRICHLET.equals(abType)) {
            t += alleleBalanceDirichletProbabilityLn(normal.hypotheses(), cancer.statistics(), i, j, alpha);
          } else if (BINOMIAL.equals(abType)) {
            t += alleleBalanceMultinomialProbabilityLn(normal.hypotheses(), cancer.statistics(), i, j, alpha);
          } else if (abType.startsWith(BAND)) {
            t += alleleBalanceBandLn(normal.hypotheses(), cancer.statistics(), i, j);
          } else if (abType.startsWith(MULTIPLE_BINOMIAL)) {
            t += alleleBalanceMultinomialPivotsLn(normal.hypotheses(), cancer.statistics(), i, j);
          }
        }
        mPosterior[i][j] = t;
      }
    }
    // After this point the special code for contamination no longer plays any part, mPosterior
    // is normal x cancer for both contaminated and uncontaminated callers.
    contraryEvidenceAdjustment(normal.statistics(), cancer.statistics());
    postConstruction();
  }

}
