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
    final int nSize = normalHyp.size(); // range is 0..nSize-1 inclusive
    final int cSize = cancerHyp.size(); // range is 1..cSize inclusive
    final double[][] q = new double[nSize][cSize + 1];
    for (int nh = 0; nh < nSize; ++nh) {
      final int a = normalCode.a(nh);
      final int b = normalCode.bc(nh);
      final int normalAlleleBits = (1 << a) | (1 << b); // for homozygous, only 1 bit gets set
      for (int ch = 1; ch <= cSize; ++ch) {
        // mu^(number of alleles different between normal hyp and cancer hyp)
        q[nh][ch] = muPowers[Integer.bitCount(normalAlleleBits ^ ch)];
      }
      q[nh] = MathUtils.renormalize(q[nh]);
    }
    return q;
  }
}
