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
