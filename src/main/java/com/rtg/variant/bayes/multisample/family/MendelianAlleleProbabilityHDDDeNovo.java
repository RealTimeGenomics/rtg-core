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
public final class MendelianAlleleProbabilityHDDDeNovo extends MendelianAlleleProbabilityDeNovo {

  /**
   * lookup table for Mendelian inheritance probability
   * indexes are mother, mother, child, child
   */
  static final double[][][][] LOOKUP = new double[4][4][6][6];
  static {
    Exam.assertEquals(4 * 4 * 6 * 6, init(LOOKUP));

    final MendelianAlleleProbabilityHDD partner = (MendelianAlleleProbabilityHDD) MendelianAlleleProbabilityHDD.SINGLETON_HD;

    // For all the parent allele combos

    // father fixed (haploid)

      for (int j = 0; j <= 1; ++j) { // mother a = 0,1
        final int mj = Math.max(0, j);

        for (int k = 0; k <= mj + 1; ++k) { // mother b = 0,1,2
        final int mk = Math.max(mj, k);

          // Now calculate de novo prob for each child hypothesis
          for (int l = 0; l <= mk + 1; ++l) { // child a = 0,1,2,3
            final int ml = Math.max(mk, l);

            for (int m = 0; m <= ml + 1; ++m) { // child a = 0,1,2,3,4

              //System.err.print("father: 0   mother: " + j + "/" + k + "  child: " + l + "/" + m + "    ");

              if (partner.getLookup(j, k, l, m) != Double.NEGATIVE_INFINITY) {  // Valid mendelian combos get de novo = 0
                //System.err.println("Mendelian, setting score to 0");
              } else {

                // De novo child hypothesis.
                double tot = 0;
                for (int n = 0; n < 5; ++n) { // child de novo a = 0,1,2,3, b = m
                  if (n != l) {
                    final double mendelian = partner.getLookup(j, k, n, m);
                    if (mendelian != Double.NEGATIVE_INFINITY) {
                      tot += Math.exp(mendelian);
                    }
                  }
                }
                for (int n = 0; n < 6; ++n) { // child de novo a = l, b = 0,1,2,3,4
                  if (n != m) {
                    final double mendelian = partner.getLookup(j, k, l, n);
                    if (mendelian != Double.NEGATIVE_INFINITY) {
                      tot += Math.exp(mendelian);
                    }
                  }
                }
                tot = tot / 2;
                //System.err.println("Non-mendelian, setting score to " + tot);
                setLookup(j, k, l, m, tot);

              }
            }
          }
        }
      }

  }

  private static int init(double[][][][] lookup) {
    final double z = Double.NEGATIVE_INFINITY;
    int cnt = 0;
    for (int i = 0; i < lookup.length; ++i) {
      for (int j = 0; j < lookup[i].length; ++j) {
        for (int k = 0; k < lookup[i][j].length; ++k) {
          for (int l = 0; l < lookup[i][j][k].length; ++l) {
            ++cnt;
            lookup[i][j][k][l] = z;
          }
        }
      }
    }
    return cnt;
  }

  private static void setLookup(int moa, int mob, int cha, int chb, double val) {
    final double log = Math.log(val);
    LOOKUP[moa][mob][cha][chb] = log;
    LOOKUP[moa][mob][chb][cha] = log;
    LOOKUP[mob][moa][cha][chb] = log;
    LOOKUP[mob][moa][chb][cha] = log;
  }

  /** Single instance for haploid-diploid parents and diploid child. */
  public static final MendelianAlleleProbability SINGLETON_HD = new MendelianAlleleProbabilityHDDDeNovo();

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
    final UniqueId uid = new UniqueId(code.rangeSize());
    final int fa0 = uid.addId(code.a(father));
    assert fa0 == 0;
    final int mo0 = uid.addId(code.a(mother));
    final int mo1 = uid.addId(code.bc(mother));
    final int ch0 = uid.addId(code.a(child));
    final int ch1 = uid.addId(code.bc(child));
    assert code.rangeSize() > 0;
    return LOOKUP[mo0][mo1][ch0][ch1] - MendelianAlleleProbabilityDiploidDeNovo.cacheLog(code.rangeSize() - 1);
  }

}
