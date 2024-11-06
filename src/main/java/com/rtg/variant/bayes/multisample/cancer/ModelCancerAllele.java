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
  private final int mDescriptionSize;

  /**
   * @param hyp used to get the probabilities in the various categories for both normal and cancer.
   * @param statistics for statistics collection
   */
  public ModelCancerAllele(final HypothesesPowerSet<D> hyp, final Statistics<?> statistics) {
    super(hyp, statistics, new NoAlleleBalance());
    mChi = (hyp.size() + 1) / 2;
    mDescriptionSize = hyp.description().size();
  }

  private ModelCancerAllele(ModelCancerAllele<D> original) {
    super(original);
    mChi = original.mChi;
    mDescriptionSize = original.mDescriptionSize;
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
    // todo I'm ignoring evidence.probability(), is that right?  -- probably, error() is the same thing but not divvied up
    // todo but less clear this is right thing to do with complex alleles ...
    final double pEr = r * evidence.pe();
    final double error = evidence.error();
    final int a = evidence.read();
    final double norm = mDescriptionSize / ((1 << mDescriptionSize) - 1 + mDescriptionSize * error);
    for (int hyp = 0; hyp < size(); ++hyp) {
      final int code = hyp + 1;
      final boolean matches = (code & (1 << a)) != 0;
      final double prob = matches ? norm / Integer.bitCount(code) : norm * error / (mChi - 1);
      if (prob <= 0.0) {
        return; // Phred scores of 0 can result in 0 probability, just skip them
      }
      // Adjust for mapQ
      final double pr = prob * rc + pEr;
      mPosteriors[hyp] = arithmetic().multiply(mPosteriors[hyp], arithmetic().prob2Poss(pr));
    }
  }

  @Override
  public HypothesisScore best(HypothesesPrior<?> hypotheses) {
    throw new UnsupportedOperationException(); // make sure we are not using this
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("Allele Cancer Model").append(LS);
    final FormatReal fmt = new FormatReal(4, 3);
    final int pad = hypotheses().maxNameLength();
    for (int hyp = 0; hyp < hypotheses().size(); ++hyp) {
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
