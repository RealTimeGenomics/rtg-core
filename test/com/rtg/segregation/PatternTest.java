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
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class PatternTest extends TestCase {

  public void test() throws MismatchingPloidyException {
    check("0/0", "0/0", "0/0", "{0/0, 0/1, 1/0, 1/1}");
    check("0/0", "0/1", "0/0", "{0/0, 1/0}");
    check("0/0", "0/1", "0/1", "{0/1, 1/1}");
    check("0/1", "0/1", "0/0", "{0/0}");
    check("0/1", "0/1", "0/1", "{0/1, 1/0}");
    check("0/1", "0/1", "1/1", "{1/1}");
    check("0/0", "1/1", "0/1", "{0/0, 0/1, 1/0, 1/1}");
    check("0/1", "2/3", "0/2", "{0/0}");
    check("0/1", "2/3", "0/3", "{0/1}");
    check("0/1", "2/3", "1/2", "{1/0}");
    check("0/1", "2/3", "1/3", "{1/1}");

    check("0/1", "1/1", "0/1", "{0/0, 0/1}");
  }

  private void check(String fa, String mo, String ch, String exp) throws MismatchingPloidyException {
    final GType fag = new GType(fa, Ploidy.DIPLOID);
    final GType mog = new GType(mo, Ploidy.DIPLOID);
    final GType chg = new GType(ch, Ploidy.DIPLOID);
    final Pattern pat = new Pattern(fag, mog, chg);
    pat.integrity();
    assertEquals(exp, pat.toString());
    //The following line generates the html table used as a comment at the beginning of Pattern
    //System.err.println(" * <tr>" + " <td>" + fa + "</td> <td>" + mo + "</td> <td>" + ch + "</td> <td>" + Utils.toBits(pat.intPattern(), 4) + "</td> <td>" + pat.intPattern() + "</td> <td>" + exp + "</td> </tr>");
  }

  public void testStringToPattern() {
    //Absolute Knowledge
    checkStringToPattern(1, false);
    checkStringToPattern(2, false);
    checkStringToPattern(4, false);
    checkStringToPattern(8, false);

    //Partial Knowledge
    checkStringToPattern(3, false);
    checkStringToPattern(5, false);
    checkStringToPattern(10, false);
    checkStringToPattern(12, false);

    //Insufficient Knowledge
    checkStringToPattern(0, true);
    checkStringToPattern(6, true);
    checkStringToPattern(7, true);
    checkStringToPattern(9, true);
    checkStringToPattern(11, true);
    checkStringToPattern(13, true);
    checkStringToPattern(14, true);
    checkStringToPattern(15, true);

    assertNull(Pattern.stringToPattern('X', '1'));
    assertNull(Pattern.stringToPattern('1', 'X'));
    assertNotNull(Pattern.stringToPattern('1', '1'));
  }

  private void checkStringToPattern(int i, boolean expNull) {
    final Pattern pat = new Pattern(i);
    final char fa = pat.faString().charAt(0);
    final char mo = pat.moString().charAt(0);
    final Pattern actual = Pattern.stringToPattern(fa, mo);
    if (expNull) {
      assertNull("" + i, actual);
    } else {
      assertNotNull("" + i, actual);
      assertEquals(i, actual.intPattern());
    }
  }

  public void testFaMoStr() {
    checkFaMo(0, "X", "X");
    checkFaMo(15, "?", "?");
    checkFaMo(1, "0", "0");
    checkFaMo(2, "0", "1");
    checkFaMo(4, "1", "0");
    checkFaMo(8, "1", "1");
    checkFaMo(3, "0", "?");
  }

  private void checkFaMo(int i, String expFa, String expMo) {
    final Pattern p = new Pattern(i);
    assertEquals(expFa, p.faString());
    assertEquals(expFa, p.sString(true));
    assertEquals(expMo, p.moString());
    assertEquals(expMo, p.sString(false));
  }

  public void testUnique() {
    final boolean[] exp = {false, true, true, false, true, false, false, false, true, false, false, false, false, false, false, false};
    for (int i = 0; i < 16; i++) {
      checkUnique(i, exp[i]);
    }
  }

  private void checkUnique(int i, boolean exp) {
    final Pattern p = new Pattern(i);
    assertEquals(exp, p.isUniqueAndValid());
  }

  public void testBit() {
    assertEquals(0, Pattern.bit(6, 0));
    assertEquals(1, Pattern.bit(6, 1));
    assertEquals(1, Pattern.bit(6, 2));
    assertEquals(0, Pattern.bit(6, 3));
  }

  public void testSet() {
    assertEquals(0, Pattern.set(0, 0));
    assertEquals(1, Pattern.set(1, 0));
    assertEquals(2, Pattern.set(1, 1));
    assertEquals(4, Pattern.set(1, 2));
    assertEquals(8, Pattern.set(1, 3));
  }

  public void testFlipInts() {
    checkFlip(0,  0, 0, 0);
    checkFlip(0, 15, 0, 0);

    checkFlip(15, 1, 0, 1);
    checkFlip(15, 1, 1, 4);
    checkFlip(15, 1, 2, 2);
    checkFlip(15, 1, 3, 8);

    checkFlip(15, 2, 0, 2);
    checkFlip(15, 2, 1, 8);
    checkFlip(15, 2, 2, 1);
    checkFlip(15, 2, 3, 4);

    checkFlip(15, 4, 0, 4);
    checkFlip(15, 4, 1, 1);
    checkFlip(15, 4, 2, 8);
    checkFlip(15, 4, 3, 2);

    checkFlip(15, 8, 0, 8);
    checkFlip(15, 8, 1, 2);
    checkFlip(15, 8, 2, 4);
    checkFlip(15, 8, 3, 1);
  }

  void checkFlip(final int a, final int b, final int f, final int exp) {
    final Pattern ap = Pattern.pattern(a);
    ap.integrity();
    final Pattern bp = Pattern.pattern(b);
    bp.integrity();
    final Pattern actual = Pattern.flipIntersect(ap, bp, f);
    if (exp == 0) {
      assertNull(actual);
    } else {
      actual.integrity();
      assertEquals(exp, actual.intPattern());
    }
  }

  public void testHash() {
    assertEquals(31, new Pattern(0).hashCode()); //regression
  }

  public void testEquals() {
    final Pattern[][] x = new Pattern[16][2];
    for (int i = 0; i < 16; i++) {
      x[i][0] = new Pattern(i);
      x[i][1] = new Pattern(i);
    }
    TestUtils.equalsHashTest(x);
  }

  public void testfaMo2Flip() {
    assertEquals(0, Pattern.faMo2Flip(false, false));
    assertEquals(1, Pattern.faMo2Flip(true, false));
    assertEquals(2, Pattern.faMo2Flip(false, true));
    assertEquals(3, Pattern.faMo2Flip(true, true));
  }

  public void testfaMo2FlipConsistent() {
    final Pattern pa = new Pattern(1);
    assertEquals("0", pa.faString());
    assertEquals("1", new Pattern(pa.flip(Pattern.faMo2Flip(true, false))).faString());
  }


}
