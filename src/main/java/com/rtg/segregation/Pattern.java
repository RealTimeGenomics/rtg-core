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

package com.rtg.segregation;

import com.rtg.reference.Ploidy;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * A single phasing pattern for one child. Each pattern is represented by 4 bits.
 * The tricky problem being dealt with is that it is unknown what the absolute phasing is and also in some cases
 * there is a constraint on the phasing of the father and mother which doesn't allow you to identify a constraint on just one of them.
 * The following table shows the alleles for the parents and child and the corresponding pattern written as a bit string, an integer and a set of allowed
 * phasings. Note that 0 and 1 have three quite distinct meanings in the alleles, the bits and the phase set.
 * <table><caption>phasing patterns</caption>
 * <tr> <th>father</th> <th>mother</th> <th>child</th> <th>bits</th> <th>int</th> <th>phase set</th> </tr>
 * <tr> <td>0_0</td> <td>0_0</td> <td>0_0</td> <td>1111</td> <td>15</td> <td>{0/0, 0/1, 1/0, 1/1}</td> </tr>
 * <tr> <td>0_0</td> <td>0_1</td> <td>0_0</td> <td>0101</td> <td>5</td> <td>{0/0, 1/0}</td> </tr>
 * <tr> <td>0_0</td> <td>0_1</td> <td>0_1</td> <td>1010</td> <td>10</td> <td>{0/1, 1/1}</td> </tr>
 * <tr> <td>0_1</td> <td>0_1</td> <td>0_0</td> <td>0001</td> <td>1</td> <td>{0/0}</td> </tr>
 * <tr> <td>0_1</td> <td>0_1</td> <td>0_1</td> <td>0110</td> <td>6</td> <td>{0/1, 1/0}</td> </tr>
 * <tr> <td>0_1</td> <td>0_1</td> <td>1_1</td> <td>1000</td> <td>8</td> <td>{1/1}</td> </tr>
 * <tr> <td>0_0</td> <td>1_1</td> <td>0_1</td> <td>1111</td> <td>15</td> <td>{0/0, 0/1, 1/0, 1/1}</td> </tr>
 * <tr> <td>0_1</td> <td>2_3</td> <td>0_2</td> <td>0001</td> <td>1</td> <td>{0/0}</td> </tr>
 * <tr> <td>0_1</td> <td>2_3</td> <td>0_3</td> <td>0010</td> <td>2</td> <td>{0/1}</td> </tr>
 * <tr> <td>0_1</td> <td>2_3</td> <td>1_2</td> <td>0100</td> <td>4</td> <td>{1/0}</td> </tr>
 * <tr> <td>0_1</td> <td>2_3</td> <td>1_3</td> <td>1000</td> <td>8</td> <td>{1/1}</td> </tr>
 * </table>
 *
 * The 0 and 1 in the phasings label groups of children with the same phasings. However these cannot be known absolutely so
 * when comparing two different phasings the four possible flips when the 0s and 1s are exchanged in the father and mother must
 * be taken into account.
 */
public class Pattern extends IntegralAbstract {

  private static final int P00 = set(1, 0);
  private static final int P01 = set(1, 1);
  private static final int P10 = set(1, 2);
  private static final int P11 = set(1, 3);

  private static final int F0 = P00 | P01;
  private static final int F1 = P10 | P11;
  private static final int M0 = P00 | P10;
  private static final int M1 = P01 | P11;

  static final int NUMBER_FLIPS = 4;
  static final int NUMBER_BITS = 4;
  private static final int SIZE_SET = 1 << NUMBER_BITS;

  private static final int[][] FLIP_TABLE = {
    {0, 1, 2, 3},
    {2, 3, 0, 1},
    {1, 0, 3, 2},
    {3, 2, 1, 0}
  };

  /**
   * Given whether father and mother have been flipped return the flip code.
   * @param fa true iff father has been flipped.
   * @param mo true iff mother has been flipped.
   * @return the flip code.
   */
  static int faMo2Flip(final boolean fa, final boolean mo) {
    return (fa ? 1 : 0) + (mo ? 2 : 0);
  }

  static int bit(final int x, final int i) {
    return (x >> i) & 1;
  }

  static int set(final int x, final int i) {
    return x << i;
  }

  /** Fast lookup cache for the transformations to the patterns by a flip. */
  static final int[][] FLIP = new int[NUMBER_FLIPS][SIZE_SET];
  static {
    for (int i = 0; i < NUMBER_FLIPS; ++i) {
      for (int j = 0; j < SIZE_SET; ++j) {
        int t = 0;
        for (int k = 0; k < NUMBER_BITS; ++k) {
          t |= set(bit(j, k), FLIP_TABLE[i][k]);
        }
        FLIP[i][j] = t;
      }
    }
  }

  /** Cache of Pattern objects given their integer values. */
  private static final Pattern[] PATTERNS = new Pattern[SIZE_SET];
  static {
    for (int i = 0; i < 16; ++i) {
      PATTERNS[i] = new Pattern(i);
    }
  }

  //Used as a utility in the tests.
  static Pattern pattern(final int i) {
    return PATTERNS[i];
  }

  /**
   * Given the output of a Pattern invert it back to the pattern. Because a "?" has lost some information
   * this can only be done for 0 and 1 values.
   * @param fa the father's phase.
   * @param mo the mother's phase.
   * @return a pattern or null if no unique pattern exists.
   */
  static Pattern stringToPattern(final char fa, final char mo) {
    //We are allowed one or the other to be ?, but not both
    if (fa == '?' && fa == mo) {
      return null;
    }
    final int moi;
    switch (mo) {
      case '0':
        moi = M0;
        break;
      case '1':
        moi = M1;
        break;
      case '?':
        moi = M0 | M1;
        break;
      case 'X':
        return null;
      default:
        throw new RuntimeException("mo=" + mo);
    }
    final int fai;
    switch (fa) {
      case '0':
        fai = F0;
        break;
      case '1':
        fai = F1;
        break;
      case '?':
        fai = F0 | F1;
        break;
      case 'X':
        return null;
      default:
        throw new RuntimeException("fa=" + fa);
    }
    return pattern(moi & fai);
  }

  /**
   * Flip b as prescribed by flip and return a pattern which results from intersecting it with
   * a. Return null if the result is empty (zero).
   * @param a first pattern to be intersected.
   * @param b second pattern to flipped and intersected.
   * @param flip specifies the flip - 0 - identity, 1 - flip the father, 2 - flip the mother, 3 - flip both father and mother.
   * @return the pattern
   */
  static Pattern flipIntersect(final Pattern a, final Pattern b, final int flip) {
    final int flipInter = a.mPattern & FLIP[flip][b.mPattern];
    if (flipInter == 0) {
      return null;
    }
    return PATTERNS[flipInter];
  }

  /**
   * Construct the integer pattern for the autosomal diploid case.
   * @param fa father.
   * @param mo mother.
   * @param a one of the alleles for a single child.
   * @param b other allele for a single child.
   * @return the pattern for a single child.
   */
  static int pat(final GType fa, final GType mo, final int a, final int b) {
    int pattern = 0;
    if (a == fa.a()) {
      if (b == mo.a()) {
        pattern |= P00;
      }
      if (b == mo.b()) {
        pattern |= P01;
      }
    }
    if (a == fa.b()) {
      if (b == mo.a()) {
        pattern |= P10;
      }
      if (b == mo.b()) {
        pattern |= P11;
      }
    }
    return pattern;
  }

  /**
   * Construct the integer pattern for X chromosome with son.
   * @param mo mother.
   * @param a only allele for son.
   * @return the pattern for son.
   */
  static int pat(final GType mo, final int a) {
    int pattern = 0;
    if (a == mo.a()) {
      pattern |= M0;
    }
    if (a == mo.b()) {
      pattern |= M1;
    }
    return pattern;
  }

  /** Stores four bits one for each of the father/mother patterns 0/0 0/1 1/0 1/1 */
  private final int mPattern;

  //For testing
  Pattern(final int pattern) {
    mPattern = pattern;
  }

  Pattern(final GType fa, final GType mo, final GType ch) {
    final int a = ch.a();
    if (fa.ploidy() == Ploidy.HAPLOID && mo.ploidy() == Ploidy.DIPLOID && ch.ploidy() == Ploidy.HAPLOID) {
      mPattern = pat(mo, a);
    } else {
      final int b = ch.b();
      mPattern = pat(fa, mo, a, b) | pat(fa, mo, b, a);
    }
  }

  int flip(final int flip) {
    return FLIP[flip][mPattern];
  }

  //Used for testing
  int intPattern() {
    return mPattern;
  }

  String sString(final boolean select) {
    return select ? faString() : moString();
  }

  String moString() {
    if (mPattern == 0) {
      return "X";
    }
    final boolean m0 = (mPattern & M0) != 0;
    final boolean m1 = (mPattern & M1) != 0;
    return m0 && m1 ? "?" : m0 ? "0" : "1";
  }

  String faString() {
    if (mPattern == 0) {
      return "X";
    }
    final boolean f0 = (mPattern & F0) != 0;
    final boolean f1 = (mPattern & F1) != 0;
    return f0 && f1 ? "?" : f0 ? "0" : "1";
  }

  boolean isUniqueAndValid() {
    if (mPattern == 0) {
      return false;
    }
    final boolean f0 = (mPattern & F0) != 0;
    final boolean f1 = (mPattern & F1) != 0;
    if (f0 && f1) {
      return false;
    }
    final boolean m0 = (mPattern & M0) != 0;
    final boolean m1 = (mPattern & M1) != 0;
    if (m0 && m1) {
      return false;
    }
    return true;
  }

  boolean isUniqueAndValidXChromosome() {
    if (mPattern == 0) {
      return false;
    }
    final boolean m0 = (mPattern & M0) != 0;
    final boolean m1 = (mPattern & M1) != 0;
    if (m0 && m1) {
      return false;
    }
    return true;
  }

  //TODO test
  @Override
  public int hashCode() {
    return mPattern + 31;
  }

  //TODO test
  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final Pattern that = (Pattern) obj;
    return this.mPattern == that.mPattern;
  }

  @Override
  public void toString(StringBuilder sb) {
    final InnerAppend is = InnerAppend.innerSet();
    sb.append(is.start());
    if ((mPattern & P00) != 0) {
      sb.append(is.inner("0/0"));
    }
    if ((mPattern & P01) != 0) {
      sb.append(is.inner("0/1"));
    }
    if ((mPattern & P10) != 0) {
      sb.append(is.inner("1/0"));
    }
    if ((mPattern & P11) != 0) {
      sb.append(is.inner("1/1"));
    }
    sb.append(is.end());
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(0 <= mPattern && mPattern <= 16);
    Exam.assertEquals(16, PATTERNS.length);
    Exam.assertEquals(4, FLIP_TABLE.length);
    Exam.assertEquals(4, FLIP.length);
    for (int i = 0; i < 4; ++i) {
      Exam.assertEquals(4, FLIP_TABLE[i].length);
      Exam.assertEquals(16, FLIP[i].length);
    }
    return true;
  }
}
