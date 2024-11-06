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
package com.rtg.variant.format;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.util.MathUtils;
import com.rtg.variant.VariantSample;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Utility methods for GL field
 */
public final class GenotypeLikelihoodUtils {
  private GenotypeLikelihoodUtils() { }

  /**
   * Produce diploid GL field values for a set of alleles
   * from VCF spec the canonical ordering is: position of hypothesis <code>j/k</code> is given by <code>F(j/k) = (k*(k+1)/2)+j</code>.
   * @param alleles the alleles to include in the likelihoods
   * @param genotypeLikelihoods the full genotype likelihoods, including hypotheses not called
   * @return array of likelihoods in the canonical order, or null if they could not be computed
   */
  static double[] diploidLikelihoods(List<String> alleles, Map<Set<String>, Double> genotypeLikelihoods) {
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final int size = alleles.size() * (alleles.size() + 1) / 2;
    final double[] likelihoods = new double[size];
    double sum = arith.zero();
    for (int j = 0; j < alleles.size(); ++j) {
      final String a = alleles.get(j);
      for (int k = j; k < alleles.size(); ++k) {
        final String b = alleles.get(k);
        // canonical ordering, see vcf spec
        final int pos = (k * (k + 1) / 2) + j;
        final Double likelihood = genotypeLikelihoods.get(VariantSample.pairSet(a, b));
        if (likelihood == null) { // Can occur if variant splitting discovers a new VA
          return null;
        }
        likelihoods[pos] = likelihood;
        sum = arith.add(sum, likelihoods[pos]);
      }
    }

    for (int i = 0; i < likelihoods.length; ++i) {
      likelihoods[i] = arith.divide(likelihoods[i], sum) / MathUtils.LOG_10;
    }
    return likelihoods;
  }

  /**
   * Produce haploid GL field values for a set of alleles
   * from VCF spec the canonical ordering is: position of hypothesis <code>j/k</code> is given by <code>F(j/k) = (k*(k+1)/2)+j</code>.
   * @param alleles the alleles to include in the likelihoods
   * @param genotypeLikelihoods the full genotype likelihoods, including hypotheses not called
   * @return array of likelihoods in the canonical order, or null if they could not be computed
   */
  static double[] haploidLikelihoods(List<String> alleles, Map<Set<String>, Double> genotypeLikelihoods) {
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final double[] likelihoods = new double[alleles.size()];
    double sum = arith.zero();
    for (int i = 0; i < alleles.size(); ++i) {
      final Double likelihood = genotypeLikelihoods.get(Collections.singleton(alleles.get(i)));
      if (likelihood == null) { // Can occur if variant splitting discovers a new VA
        return null;
      }
      likelihoods[i] = likelihood;
      sum = arith.add(sum, likelihoods[i]);
    }
    for (int i = 0; i < alleles.size(); ++i) {
      likelihoods[i] = arith.divide(likelihoods[i], sum) / MathUtils.LOG_10;
    }
    return likelihoods;
  }
}
