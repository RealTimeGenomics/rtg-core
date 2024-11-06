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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.variant.bayes.Code;

/**
 * Calculates <code>P'dagger'</code> for the case where one parent chromosomes is diploid and the other haploid
 * - X chromosome for mammalian genomes.
 * Handles arbitrary hypothesis sets (eg for complex calling).
 * See the <code>multiScoring.tex</code> documentation for Consistent Mendelian pedigrees
 */
//this version has been modified to handle arbitrary hypothesis sets (eg for complex calling).
public class AlleleSetProbabilityHaploid implements AlleleProbability {

  static final double[] G = {1, 14, 36, 24};

  static double gSum() {
    return G[0] * 4 + G[1] * 6 + G[2] * 4 + G[3];
  }

  /**
   * lookup table for the probability of selecting the specified haploid sets from
   * the alleles present in the set + ref
   * indexes are reference, mother, mother
   */
  static final double[][][] LOOKUP = new double[4][2][3];
  static {
    for (byte i = 0; i < 4; ++i) {
      for (int k = 0; k < 2; ++k) {
        for (int l = 0; l < 3; ++l) {
          setLookup(i, k, l);
        }
      }
    }
  }

  private static void setLookup(int r, int m0, int m1) {
    final UniqueId uid = new UniqueId(5);
    final int hm = m0 == m1 ? 1 : 2;

    uid.addId(0);
    uid.addId(m0);
    uid.addId(m1);
    final int n = uid.numberIdsSoFar();
    final int rid = uid.addId(r);
    final int f = rid == n ? 1 : n;
    final int nn = uid.numberIdsSoFar();
    final double v = G[nn - 1] / (hm * f);
    //if (r == 0 && f0 == 0 && f1 == 0 && m0 == 1 && m1 == 2) {
    //System.err.println("v=" + v + " n=" + n + " f=" + f + " hf=" + hf + " hm=" + hm + " rid=" + rid);
    //}
    LOOKUP[r][m0][m1] = -Math.log(v);
  }

  /** One instance that handles the haploid-diploid case. */
  public static final AlleleProbability SINGLETON_HD = new AlleleSetProbabilityHaploid();

  /** One instance that handles the haploid-diploid case. */
  public static final AlleleProbability SINGLETON_DH = new AlleleSetProbabilityHaploid() {
    @Override
    public double probabilityLn(Code code, int ref, int father, int mother) {
      return SINGLETON_HD.probabilityLn(code, ref, mother, father);
    }
  };

  @Override
  public double probabilityLn(Code code, int ref, int father, int mother) {
    assert code.homozygous(father);
    final UniqueId uid = new UniqueId(code.rangeSize());
    final int fa = uid.addId(code.a(father));
    assert fa == 0;
    final int mo0 = uid.addId(code.a(mother));
    final int mo1 = uid.addId(code.bc(mother));
    final int r = ref == -1 ? 0 : uid.addId(ref); // If reference is unknown, we cannot treat it as another allele in the set.
    return LOOKUP[r][mo0][mo1];
  }

}
