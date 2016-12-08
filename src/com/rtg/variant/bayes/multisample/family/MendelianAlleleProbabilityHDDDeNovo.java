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
