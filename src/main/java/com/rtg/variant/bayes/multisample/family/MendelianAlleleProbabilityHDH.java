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

import com.rtg.util.integrity.Exam;
import com.rtg.variant.bayes.Code;

/**
 * Calculates the probability of a particular child allele combination given parents.
 * Parents are one diploid, one haploid and child is haploid.
 */
public final class MendelianAlleleProbabilityHDH extends MendelianAlleleProbabilityNormal {

  /**
   * lookup table for Mendelian inheritance probability
   * indexes are mother, child
   */
  static final double[][] LOOKUP = new double[2][];
  static {

    LOOKUP[0] = new double[1];

    LOOKUP[1] = new double[2];

    Exam.assertEquals(3, init(LOOKUP));

    setLookup(0, 0, 1.0);

    setLookup(1, 0, 0.5);
    setLookup(1, 1, 0.5);
  }

  private static int init(double[][] lookup) {
    final double z = Double.NEGATIVE_INFINITY;
    int cnt = 0;
    for (int i = 0; i < lookup.length; ++i) {
      for (int j = 0; j < lookup[i].length; ++j) {
        ++cnt;
        lookup[i][j] = z;
      }
    }
    return cnt;
  }

  private static void setLookup(int mob, int cha, double val) {
    final double log = Math.log(val);
    LOOKUP[mob][cha] = log;
  }

  double getLookup(int mob, int cha) {
    final double[] parents = LOOKUP[mob];
    if (parents.length > cha) {
      return parents[cha];
    } else {
      return Double.NEGATIVE_INFINITY;
    }
  }


  /** Single instance for haploid-diploid parents and haploid child. */
  public static final MendelianAlleleProbability SINGLETON_HD = new MendelianAlleleProbabilityHDH();

  /** Single instance for diploid-haploid parents and haploid child. */
  public static final MendelianAlleleProbability SINGLETON_DH = new MendelianAlleleProbabilityNormal() {
    @Override
    public double probabilityLn(Code code, int father, int mother, int child) {
      return SINGLETON_HD.probabilityLn(code, mother, father, child);
    }
  };

  @Override
  public double probabilityLn(Code code, int father, int mother, int child) {
    assert code.homozygous(father);
    assert code.homozygous(child);
    //System.err.println("father=" + father + " mother=" + mother + " child=" + child + " code=" + code);
    final UniqueId uid = new UniqueId(code.rangeSize());
    final int mo0 = uid.addId(code.a(mother));
    assert mo0 == 0;
    final int mo1 = uid.addId(code.bc(mother));
    final int ch0 = uid.id(code.a(child));
    if (ch0 < 0) {
      return Double.NEGATIVE_INFINITY;
    }

    return LOOKUP[mo1][ch0];
  }

}
