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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.variant.bayes.Code;

/**
 * Calculates <code>P'dagger'</code> for the case where both parent chromosomes are diploid - the normal autosomal
 * case, including arbitrary hypothesis sets (eg for complex calling).
 * See the <code>multiScoring.tex</code> documentation for Consistent Mendelian pedigrees
 */
public final class AlleleSetProbabilityDiploid implements AlleleProbability {
  /** Utility class constructor */
  private AlleleSetProbabilityDiploid() { }

  private static final double[] G = {1, 30, 150, 240, 120};

  static double gSum() {
    return G[0] * 4 + G[1] * 6 + G[2] * 4 + G[3];
  }

  /**
   * lookup table for the probability of selecting the specified haploid sets from
   * the alleles present in the set + ref
   * indexes are reference, father-b, mother-a, mother-a
   */
  static final double[][][][] LOOKUP = new double[5][2][3][4];
  static {
    for (byte i = 0; i < 5; ++i) {
      for (int k = 0; k < 2; ++k) {
        for (int l = 0; l < 3; ++l) {
          for (int m = 0; m < 4; ++m) {
            setLookup(i, k, l, m);
          }
        }
      }
    }
  }

  private static void setLookup(int r, int f1, int m0, int m1) {
    final UniqueId uid = new UniqueId(5);
    final int hf = 0 == f1 ? 1 : 2;
    final int hm = m0 == m1 ? 1 : 2;

    uid.addId(0);
    uid.addId(f1);
    uid.addId(m0);
    uid.addId(m1);
    final int n = uid.numberIdsSoFar();
    final int rid = uid.addId(r);
    final int f = rid == n ? 1 : n;
    final int nn = uid.numberIdsSoFar();
    final double v = G[nn - 1] / (hf * hm * f);
    //if (r == 0 && f0 == 0 && f1 == 0 && m0 == 1 && m1 == 2) {
    //System.err.println("v=" + v + " n=" + n + " f=" + f + " hf=" + hf + " hm=" + hm + " rid=" + rid);
    //}
    LOOKUP[r][f1][m0][m1] = -Math.log(v);
  }

  /** One instance that handles the diploid-diploid case. */
  public static final AlleleProbability SINGLETON = new AlleleSetProbabilityDiploid();

  @Override
  public double probabilityLn(Code code, int ref, int father, int mother) {
    final UniqueId uid = new UniqueId(code.rangeSize());
    final int fa0 = uid.addId(code.a(father));
    assert fa0 == 0;
    final int fa1 = uid.addId(code.bc(father));
    final int mo0 = uid.addId(code.a(mother));
    final int mo1 = uid.addId(code.bc(mother));
    final int r = ref == -1 ? 0 : uid.addId(ref); // If reference is unknown, we cannot treat it as another allele in the set.
    return LOOKUP[r][fa1][mo0][mo1];
  }

}
