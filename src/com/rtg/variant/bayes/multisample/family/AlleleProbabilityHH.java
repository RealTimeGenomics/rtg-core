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
 */
public final class AlleleProbabilityHH  implements AlleleProbability {
  /** Utility class constructor */
  private AlleleProbabilityHH() { }

  private static final double[] G = {1, 6, 6};

  static double gSum() {
    return G[0] * 4 + G[1] * 6 + G[2] * 4;
  }

  /**
   * lookup table for the probability of selecting the specified haploid sets from
   * the alleles present in the set + ref
   * indexes are mother-a, reference
   */
  static final double[][] LOOKUP = new double[2][3];
  static {
    for (byte i = 0; i < 2; ++i) {
      for (int k = 0; k < 3; ++k) {
        setLookup(i, k);
      }
    }
  }

  private static void setLookup(int m0, int r) {
    final UniqueId uid = new UniqueId(5);
    uid.addId(0);
    uid.addId(m0);
    final int n = uid.numberIdsSoFar(); // |H_f U H_m|
    uid.addId(r);
    final int nn = uid.numberIdsSoFar(); // |A|
    final int f = nn == n ? 1 : n;
    final double v = G[nn - 1] / f;
    //if (r == 0 && f0 == 0 && f1 == 0 && m0 == 1 && m1 == 2) {
    //System.err.println("v=" + v + " n=" + n + " f=" + f + " hf=" + hf + " hm=" + hm + " rid=" + rid);
    //}
    LOOKUP[m0][r] = -Math.log(v);
  }

  /** One instance that handles the haploid-haploid case. */
  public static final AlleleProbability SINGLETON = new AlleleProbabilityHH();

  @Override
  public double probabilityLn(Code code, int ref, int father, int mother) {
    assert code.homozygous(father);
    assert code.homozygous(mother);
    final UniqueId uid = new UniqueId(code.rangeSize());
    final int fa0 = uid.addId(code.a(father));
    assert fa0 == 0;
    final int mo0 = uid.addId(code.a(mother));
    final int r = ref == -1 ? 0 : uid.addId(ref); // If reference is no hypothesis, we cannot treat it as another allele in the set.
    return LOOKUP[mo0][r];
  }
}
