/*
 * Copyright (c) 2017. Real Time Genomics Limited.
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
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.HypothesesPowerSet;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.Statistics;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.snp.HypothesesPrior;

/**
 * Accumulates posteriors for the cancer side of a cancer sample.
 * @param <D> description type
 */
public class ModelCancerAllele<D extends Description> extends Model<D> {

  /** Number of hypotheses containing any specific allele. */
  private final int mChi;

  /**
   * @param hyp used to get the probabilities in the various categories for both normal and cancer.
   * @param statistics for statistics collection
   */
  public ModelCancerAllele(final HypothesesPowerSet<D> hyp, final Statistics<?> statistics) {
    super(hyp, statistics, new NoAlleleBalance());
    mChi = (hyp.size() + 1) / 2;
  }

  private ModelCancerAllele(ModelCancerAllele<D> original) {
    super(original);
    mChi = original.mChi;
  }

  @Override
  public void freeze() {
    // Avoid doing the allele balance calculation.
    assert !mFrozen : "Should only freeze once";
    mFrozen = true;
  }

  @Override
  public void increment(final EvidenceInterface evidence) {
    incrementStatistics(evidence);
    if (ambiguityShortCircuit(evidence)) {
      return; // todo do we want this?
    }
    final double r = evidence.mapError();
    final double rc = 1.0 - r;
    if (rc <= 0.0) {
      return; // Avoid case where mapq is 0 which gives a NaN
    }
    // todo I'm ignoring evidence.probability(), is that right?
    // todo but less clear this is right thing to do with complex alleles ...
    final double pEr = r * evidence.pe();
    final double error = evidence.error();
    final int a = evidence.read();
    for (int hyp = 0; hyp < size(); ++hyp) {
      final boolean matches = ((hyp + 1) & (1 << a)) != 0;
      final double prob = matches ? (1 - error) / mChi : error / (mChi - 1);
      if (prob <= 0.0) {
        return; // Phred scores of 0 can result in 0 probabilty, just skip them
      }
      // Adjust for mapQ
      final double pr = prob * rc + pEr;
      mPosteriors[hyp] = arithmetic().multiply(mPosteriors[hyp], arithmetic().prob2Poss(pr));
    }
  }

  // xxx todo probably need to override best() since hypothesis number is non-standard


  @Override
  public HypothesisScore best(HypothesesPrior<?> hypotheses) {
    throw new UnsupportedOperationException(); // make sure we are not using this
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Allele Cancer Model");
    final FormatReal fmt = new FormatReal(4, 3);
    final int pad = hypotheses().maxNameLength();
    for (int hyp = 1; hyp <= hypotheses().size(); ++hyp) {
      sb.append(StringUtils.padLeft(hypotheses().name(hyp), pad));
      sb.append(fmt.format(arithmetic().poss2Ln(mPosteriors[hyp])));
      sb.append(LS);
    }
  }

  @Override
  public ModelCancerAllele<D> copy() {
    return new ModelCancerAllele<>(this);
  }
}
