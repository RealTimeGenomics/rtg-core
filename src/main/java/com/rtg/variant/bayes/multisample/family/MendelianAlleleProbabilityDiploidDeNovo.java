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

import com.rtg.util.MathUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.variant.bayes.Code;

/**
 * Calculates the probability of a particular child allele combination given parents
 */
public final class MendelianAlleleProbabilityDiploidDeNovo extends MendelianAlleleProbabilityDeNovo {

  static final int CACHE_SIZE = 101;
  static final double[] CACHED_LOGS = new double[CACHE_SIZE];
  static {
    for (int i = 1; i < CACHE_SIZE; ++i) {
      CACHED_LOGS[i] = Math.log(i);
    }
  }

  static double cacheLog(int val) {
    if (val < CACHE_SIZE) {
      return CACHED_LOGS[val];
    }
    return MathUtils.log(val);
  }

  /** Utility class private constructor */
  private MendelianAlleleProbabilityDiploidDeNovo() { }

  /**
   * lookup table for Mendelian inheritance probability
   * indexes are father, mother, mother, child, child
   */
  static final double[][][][][] LOOKUP = new double[2][4][4][6][6];
  static {
    Exam.assertEquals(2 * 4 * 4 * 6 * 6, init(LOOKUP));

    final MendelianAlleleProbabilityDiploid partner = (MendelianAlleleProbabilityDiploid) MendelianAlleleProbabilityDiploid.SINGLETON;

    // For all the parent allele combos
    for (int i = 0; i < 2; ++i) { // father a = 0, father b = i

      for (int j = 0; j <= i + 1; ++j) { // existing allele count + 1, mother a = 0,1,2
        final int mj = Math.max(i, j); // Track highest allele

        for (int k = 0; k <= mj + 1; ++k) { // mother b = 0,1,2,3
        final int mk = Math.max(mj, k);

          // Now calculate de novo prob for each child hypothesis
          for (int l = 0; l <= mk + 1; ++l) { // child a = 0,1,2,3,4
            final int ml = Math.max(mk, l);

            for (int m = 0; m <= ml + 1; ++m) { // child a = 0,1,2,3,4,5

              //System.err.print("father: 0/" + i + "  mother: " + j + "/" + k + "  child: " + l + "/" + m + "    ");

              if (partner.getLookup(i, j, k, l, m) != Double.NEGATIVE_INFINITY) {  // Valid mendelian combos get de novo = 0
                //System.err.println("Mendelian, setting score to 0");
              } else {

                // De novo child hypothesis.
                double tot = 0;
                for (int n = 0; n < 5; ++n) { // child de novo a = 0,1,2,3,4, b = m
                  if (n != l) {
                    final double mendelian = partner.getLookup(i, j, k, n, m);
                    if (mendelian != Double.NEGATIVE_INFINITY) {
                      tot += Math.exp(mendelian);
                    }
                  }
                }
                for (int n = 0; n < 6; ++n) { // child de novo a = l, b = 0,1,2,3,4,5
                  if (n != m) {
                    final double mendelian = partner.getLookup(i, j, k, l, n);
                    if (mendelian != Double.NEGATIVE_INFINITY) {
                      tot += Math.exp(mendelian);
                    }
                  }
                }
                tot = tot / 2;
                //System.err.println("Non-mendelian, setting score to " + tot);
                setLookup(i, j, k, l, m, tot);

              }
            }
          }
        }
      }
    }
  }

  private static int init(double[][][][][] array) {
    final double z = Double.NEGATIVE_INFINITY;
    int cnt = 0;
    for (int i = 0; i < array.length; ++i) {
      for (int j = 0; j < array[i].length; ++j) {
        for (int k = 0; k < array[i][j].length; ++k) {
          for (int l = 0; l < array[i][j][k].length; ++l) {
            for (int m = 0; m < array[i][j][k][l].length; ++m) {
              ++cnt;
              array[i][j][k][l][m] = z;
            }
          }
        }
      }
    }
    return cnt;
  }

  private static void setLookup(int fab, int moa, int mob, int cha, int chb, double val) {
    final double log = Math.log(val);
    LOOKUP[fab][moa][mob][cha][chb] = log;
    LOOKUP[fab][moa][mob][chb][cha] = log;
    LOOKUP[fab][mob][moa][cha][chb] = log;
    LOOKUP[fab][mob][moa][chb][cha] = log;
  }

  /** Single instance for diploid-diploid case. */
  public static final MendelianAlleleProbability SINGLETON = new MendelianAlleleProbabilityDiploidDeNovo();

  @Override
  public double probabilityLn(Code code, int father, int mother, int child) {
    assert code.rangeSize() > 0;
    final UniqueId uid = new UniqueId(code.rangeSize());
    final int fa0 = uid.addId(code.a(father));
    assert fa0 == 0;
    final int fa1 = uid.addId(code.bc(father));
    final int mo0 = uid.addId(code.a(mother));
    final int mo1 = uid.addId(code.bc(mother));
    final int ch0 = uid.addId(code.a(child));
    final int ch1 = uid.addId(code.bc(child));
    return LOOKUP[fa1][mo0][mo1][ch0][ch1] - MendelianAlleleProbabilityDiploidDeNovo.cacheLog(code.rangeSize() - 1);
  }

}
