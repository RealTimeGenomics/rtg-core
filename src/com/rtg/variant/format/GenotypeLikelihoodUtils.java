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
package com.rtg.variant.format;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.reference.Ploidy;
import com.rtg.util.MathUtils;
import com.rtg.variant.VariantSample;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Utility methods for GL field
 */
public final class GenotypeLikelihoodUtils {
  private GenotypeLikelihoodUtils() { }
  private static double[] diploidLikelihoods(List<String> calls, Map<Set<String>, Double> genotypeLikelihoods) {
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final int size = calls.size() * (calls.size() + 1) / 2;
    final double[] likelihoods = new double[size];
    double sum = arith.zero();
    for (int j = 0; j < calls.size(); j++) {
      final String a = calls.get(j);
      for (int k = j; k < calls.size(); k++) {
        final String b = calls.get(k);
        // canonical ordering, see vcf spec
        final int pos = (k * (k + 1) / 2) + j;
        final Double likelihood = genotypeLikelihoods.get(VariantSample.pairSet(a, b));
        if (likelihood == null) {
          return null; // XXX we may not have computed the likelihood for the VA allele
        }
        likelihoods[pos] = likelihood;
        sum = arith.add(sum, likelihoods[pos]);
      }
    }

    for (int i = 0; i < likelihoods.length; i++) {
      likelihoods[i] = arith.divide(likelihoods[i], sum) / MathUtils.LOG_10;
    }
    return likelihoods;
  }

  private static double[] haploidLikelihoods(List<String> calls, Map<Set<String>, Double> genotypeLikelihoods) {
    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final double[] likelihoods = new double[calls.size()];
    double sum = arith.zero();
    for (int i = 0; i < calls.size(); i++) {
      final Double likelihood = genotypeLikelihoods.get(Collections.singleton(calls.get(i)));
      if (likelihood == null) {
        return null; // XXX we may not have computed the likelihood for the VA allele
      }
      likelihoods[i] = likelihood;
      sum = arith.add(sum, likelihoods[i]);
    }
    for (int i = 0; i < calls.size(); i++) {
      likelihoods[i] = arith.divide(likelihoods[i], sum) / MathUtils.LOG_10;
    }
    return likelihoods;
  }

  /**
   * Produce GL field values for a set of alleles
   * from VCF spec the canonical ordering is: position of hypothesis <code>j/k</code> is given by <code>F(j/k) = (k*(k+1)/2)+j</code>.
   * @param sample the sample to create likelihoods for
   * @param alleles the alleles to include in the likelihoods
   * @return array of likelihoods in the canonical order
   */
  public static double[] getLikelihoods(VariantSample sample, List<String> alleles) {
    final double[] likelihoods;
    final Map<Set<String>, Double> genotypeLikelihoods = sample.getGenotypeLikelihoods();
    if (genotypeLikelihoods == null || genotypeLikelihoods.size() == 0) {
      return null;
    }
    if (sample.getPloidy() == Ploidy.HAPLOID) {
      likelihoods = haploidLikelihoods(alleles, genotypeLikelihoods);
    } else if (sample.getPloidy() == Ploidy.DIPLOID) {
      likelihoods = diploidLikelihoods(alleles, genotypeLikelihoods);
    } else {
      likelihoods = null;
    }
    return likelihoods;
  }

}
