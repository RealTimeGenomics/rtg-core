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
 * Calculates the probability of a particular child allele combination given parents
 */
public final class MendelianAlleleProbabilityHDHDeNovo extends MendelianAlleleProbabilityDeNovo {

  /**
   * lookup table for Mendelian inheritance probability
   * indexes are mother, child
   */
  static final double[][] LOOKUP = new double[2][3];
  static {
    Exam.assertEquals(2 * 3, init(LOOKUP));

    final MendelianAlleleProbabilityHDH partner = (MendelianAlleleProbabilityHDH) MendelianAlleleProbabilityHDH.SINGLETON_HD;

    // For all the parent allele combos

    // father fixed (haploid)

    // mother fixed (0)

    for (int k = 0; k <= 1; ++k) { // mother b = 0,1

      // Now calculate de novo prob for each child hypothesis
      for (int l = 0; l <= k + 1; ++l) { // child a = 0,1,2

        //System.err.print("father: 0   mother: 0/" + k + "  child: " + l + "    ");

        if (partner.getLookup(k, l) != Double.NEGATIVE_INFINITY) {  // Valid mendelian combos get de novo = 0
          //System.err.println("Mendelian, setting score to 0");
        } else {

          // De novo child hypothesis.
          double tot = 0;
          for (int n = 0; n < 3; ++n) { // child de novo a = 0,1,2
            if (n != l) {
              final double mendelian = partner.getLookup(k, n);
              if (mendelian != Double.NEGATIVE_INFINITY) {
                tot += Math.exp(mendelian);
              }
            }
          }
          //System.err.println("Non-mendelian, setting score to " + tot);
          setLookup(k, l, tot);

        }
      }
    }
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

  /** Single instance for haploid-diploid parents and diploid child. */
  public static final MendelianAlleleProbability SINGLETON_HD = new MendelianAlleleProbabilityHDHDeNovo();

  /** Single instance for diploid-haploid parents and diploid child. */
  public static final MendelianAlleleProbability SINGLETON_DH = new MendelianAlleleProbabilityDeNovo() {
    @Override
    public double probabilityLn(Code code, int father, int mother, int child) {
      return SINGLETON_HD.probabilityLn(code, mother, father, child);
    }
  };


  @Override
  public double probabilityLn(Code code, int father, int mother, int child) {
    assert code.homozygous(father);
    assert code.homozygous(child);
    final UniqueId uid = new UniqueId(code.rangeSize());
    final int mo0 = uid.addId(code.a(mother));
    assert mo0 == 0;
    final int mo1 = uid.addId(code.bc(mother));
    final int ch0 = uid.addId(code.a(child));
    assert code.rangeSize() > 0;
    return LOOKUP[mo1][ch0] - MendelianAlleleProbabilityDiploidDeNovo.cacheLog(code.rangeSize() - 1);
  }

}
