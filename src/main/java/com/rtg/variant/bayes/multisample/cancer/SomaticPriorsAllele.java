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

import com.rtg.util.MathUtils;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.HypothesesPowerSet;

/**
 * Q matrix calculation for the allele based cancer caller
 */
final class SomaticPriorsAllele {

  private SomaticPriorsAllele() { }

  /**
   * Compute a Q matrix for the allele based cancer caller.  Essentially the number of allele
   * differences between the normal and cancer hypotheses.  This will work for both haploid
   * and diploid normal hypotheses.
   * @param mu probability of a somatic mutation.
   * @param normalHyp normal hypotheses
   * @param cancerHyp cancer allele based hypotheses
   * @param <D> description type
   * @return Q transition matrix
   */
  static <D extends Description> double[][] makeQ(final double mu, final Hypotheses<D> normalHyp, final HypothesesPowerSet<D> cancerHyp) {
    // Details of this calculation depends on the ordering of alleles in HypothesesPowerSet

    // Essentially each difference between the normal an cancer alleles gets a penalty
    // of mu.  This may need tuning, for example P(Y|X) = mu^2, but P(XY|X) = mu.

    final int descSize = normalHyp.description().size();
    // Precompute powers of mu
    final double[] muPowers = new double[descSize + 1];
    muPowers[0] = 1;
    for (int k = 1; k < muPowers.length; ++k) {
      muPowers[k] = muPowers[k - 1] * mu;
    }
    // Matrix orientation is q[normal][cancer]
    final Code normalCode = normalHyp.code();
    final int nSize = normalHyp.size();
    final int cSize = cancerHyp.size();
    final double[][] q = new double[nSize][cSize];
    for (int nh = 0; nh < nSize; ++nh) {
      final int a = normalCode.a(nh);
      final int b = normalCode.bc(nh);
      final int normalAlleleBits = (1 << a) | (1 << b); // for homozygous, only 1 bit gets set
      for (int ch = 0; ch < cSize; ++ch) {
        // mu^(number of alleles different between normal hyp and cancer hyp)
        q[nh][ch] = muPowers[Integer.bitCount(normalAlleleBits ^ (ch + 1))]; // ch + 1 to get bit coding
      }
      q[nh] = MathUtils.renormalize(q[nh]);
    }
    return q;
  }
}
