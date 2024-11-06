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
