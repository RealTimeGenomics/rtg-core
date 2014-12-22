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
package com.rtg.variant.cnv.region;

import static com.rtg.util.StringUtils.LS;

import java.util.SortedSet;
import java.util.TreeSet;

import com.rtg.util.integrity.Exam.ExamException;

import junit.framework.TestCase;

/**
 */
public class ComplexCnvRegionTest extends TestCase {

  public void testCheck() {
    final SortedSet<AbstractCnvRegion> set = new TreeSet<>();
    set.add(new SimpleCnvRegion(3, 3));
    set.add(new SimpleCnvRegion(7, 9));
    final ComplexCnvRegion re = new ComplexCnvRegion(3, 12, set);
    re.checkSubRegion(new SimpleCnvRegion(2, 3));
    re.checkSubRegion(new SimpleCnvRegion(5, 10));
    re.checkSubRegion(new SimpleCnvRegion(2, 13));
    try {
      re.checkSubRegion(new SimpleCnvRegion(1, 2));
      fail();
    } catch (final ExamException e) {
      // expected
    }
    try {
      re.checkSubRegion(new SimpleCnvRegion(13, 14));
      fail();
    } catch (final ExamException e) {
      // expected
    }
  }

  public void testBad0() {
    final SortedSet<AbstractCnvRegion> set = new TreeSet<>();
    try {
      new ComplexCnvRegion(0, 12, set);
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }
  }

  public void testBad1() {
    final SortedSet<AbstractCnvRegion> set = new TreeSet<>();
    set.add(new SimpleCnvRegion(0, 0));
    try {
      new ComplexCnvRegion(0, 12, set);
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }
  }

  public void test2() {
    final SortedSet<AbstractCnvRegion> set = new TreeSet<>();
    set.add(new SimpleCnvRegion(0, 1));
    set.add(new SimpleCnvRegion(7, 9));
    final ComplexCnvRegion re = new ComplexCnvRegion(0, 12, set);
    re.globalIntegrity();
    assertFalse(re.isInRegion(Integer.MIN_VALUE));
    assertFalse(re.isInRegion(-1));
    assertTrue(re.isInRegion(0));
    assertFalse(re.isInRegion(1));
    assertFalse(re.isInRegion(2));
    assertFalse(re.isInRegion(3));
    assertFalse(re.isInRegion(4));
    assertFalse(re.isInRegion(5));
    assertFalse(re.isInRegion(6));
    assertTrue(re.isInRegion(7));
    assertTrue(re.isInRegion(8));
    assertFalse(re.isInRegion(9));
    assertFalse(re.isInRegion(10));
    assertFalse(re.isInRegion(11));
    assertFalse(re.isInRegion(12));
    assertFalse(re.isInRegion(13));
    assertFalse(re.isInRegion(Integer.MAX_VALUE));

    assertEquals(0, re.getStart());
    assertEquals(12, re.getEnd());

    for (int i = 0; i <= 12; i++) {
      checkBox(re, i);
    }
    final String exp = ""
      + "ComplexRegion start=0 end=12 length=4" + LS
      + "First SimpleRegion start=0 end=1" + LS
      + "Last SimpleRegion start=7 end=9" + LS
      ;
    assertEquals(exp, re.toString());
  }

  public void test3() {
    final SortedSet<AbstractCnvRegion> set = new TreeSet<>();
    set.add(new SimpleCnvRegion(0, 1));
    set.add(new SimpleCnvRegion(3, 5));
    set.add(new SimpleCnvRegion(7, 9));
    final ComplexCnvRegion re = new ComplexCnvRegion(0, 12, set);
    re.globalIntegrity();
    assertFalse(re.isInRegion(Integer.MIN_VALUE));
    assertFalse(re.isInRegion(-1));
    assertTrue(re.isInRegion(0));
    assertFalse(re.isInRegion(1));
    assertFalse(re.isInRegion(2));
    assertTrue(re.isInRegion(3));
    assertTrue(re.isInRegion(4));
    assertFalse(re.isInRegion(5));
    assertFalse(re.isInRegion(6));
    assertTrue(re.isInRegion(7));
    assertTrue(re.isInRegion(8));
    assertFalse(re.isInRegion(9));
    assertFalse(re.isInRegion(10));
    assertFalse(re.isInRegion(11));
    assertFalse(re.isInRegion(12));
    assertFalse(re.isInRegion(13));
    assertFalse(re.isInRegion(Integer.MAX_VALUE));

    assertEquals(0, re.getStart());
    assertEquals(12, re.getEnd());

    for (int i = 0; i <= 12; i++) {
      checkBox(re, i);
    }
    final String exp = ""
      + "ComplexRegion start=0 end=12 length=6" + LS
      + "First SimpleRegion start=0 end=1" + LS
      + "[1]" + LS
      + "SimpleRegion start=3 end=5" + LS
      + "[2]" + LS
      + "SimpleRegion start=3 end=5" + LS
      + "Last SimpleRegion start=7 end=9" + LS
      ;
    assertEquals(exp, re.toString());
  }

  public void test6() {
    final SortedSet<AbstractCnvRegion> set = new TreeSet<>();
    set.add(new SimpleCnvRegion(0, 2));
    set.add(new SimpleCnvRegion(3, 5));
    set.add(new SimpleCnvRegion(3, 6));
    set.add(new SimpleCnvRegion(3, 7));
    set.add(new SimpleCnvRegion(7, 9));
    set.add(new SimpleCnvRegion(11, 15));
    final ComplexCnvRegion re = new ComplexCnvRegion(1, 12, set);
    re.globalIntegrity();
    assertFalse(re.isInRegion(Integer.MIN_VALUE));
    assertFalse(re.isInRegion(-5));
    assertFalse(re.isInRegion(-4));
    assertFalse(re.isInRegion(-3));
    assertFalse(re.isInRegion(-2));
    assertFalse(re.isInRegion(-1));
    assertFalse(re.isInRegion(0));
    assertTrue(re.isInRegion(1));
    assertFalse(re.isInRegion(2));
    assertTrue(re.isInRegion(3));
    assertTrue(re.isInRegion(4));
    assertTrue(re.isInRegion(5));
    assertTrue(re.isInRegion(6));
    assertTrue(re.isInRegion(7));
    assertTrue(re.isInRegion(8));
    assertFalse(re.isInRegion(9));
    assertFalse(re.isInRegion(10));
    assertTrue(re.isInRegion(11));
    assertTrue(re.isInRegion(12));
    assertFalse(re.isInRegion(13));
    assertFalse(re.isInRegion(14));
    assertFalse(re.isInRegion(15));
    assertFalse(re.isInRegion(16));
    assertFalse(re.isInRegion(17));
    assertFalse(re.isInRegion(18));
    assertFalse(re.isInRegion(Integer.MAX_VALUE));

    assertEquals(1, re.getStart());
    assertEquals(12, re.getEnd());

    for (int i = 2; i <= 12; i++) {
      checkBox(re, i);
    }
    final String exp = ""
      + "ComplexRegion start=1 end=12 length=12" + LS
      + "First SimpleRegion start=0 end=2" + LS
      + "[2]" + LS
      + "ComplexRegion start=3 end=3 length=6" + LS
      + "First SimpleRegion start=3 end=5" + LS
      + "[0]" + LS
      + "SimpleRegion start=3 end=6" + LS
      + "Last SimpleRegion start=3 end=7" + LS
      + "" + LS
      + "[3]" + LS
      + "ComplexRegion start=4 end=4 length=6" + LS
      + "First SimpleRegion start=3 end=5" + LS
      + "[0]" + LS
      + "SimpleRegion start=3 end=6" + LS
      + "Last SimpleRegion start=3 end=7" + LS
      + "" + LS
      + "[4]" + LS
      + "ComplexRegion start=5 end=5 length=6" + LS
      + "First SimpleRegion start=3 end=5" + LS
      + "[0]" + LS
      + "SimpleRegion start=3 end=6" + LS
      + "Last SimpleRegion start=3 end=7" + LS
      + "" + LS
      + "[5]" + LS
      + "ComplexRegion start=6 end=6 length=4" + LS
      + "First SimpleRegion start=3 end=6" + LS
      + "Last SimpleRegion start=3 end=7" + LS
      + "" + LS
      + "[6]" + LS
      + "ComplexRegion start=7 end=7 length=4" + LS
      + "First SimpleRegion start=3 end=7" + LS
      + "Last SimpleRegion start=7 end=9" + LS
      + "" + LS
      + "[7]" + LS
      + "SimpleRegion start=7 end=9" + LS
      + "[8]" + LS
      + "SimpleRegion start=7 end=9" + LS
      + "Last SimpleRegion start=11 end=15" + LS
      ;
    assertEquals(exp, re.toString());
  }

  private void checkBox(final ComplexCnvRegion r, final int index) {
    final int box = r.box(index);
    final int x = r.index(box);
    assertTrue("index=" + index + " box=" + box + " x=" + x, x <= index);
    final int b = r.box(x);
    assertEquals("index=" + index + " box=" + box + " x=" + x + " b=" + b, box, b);
  }

  public void testOverflow() {
    final SortedSet<AbstractCnvRegion> set = new TreeSet<>();
    for (int i = 0; i < 2640; i++) {
      set.add(new SimpleCnvRegion(i + 90, i + 91));
    }
    set.add(new SimpleCnvRegion(804577, 804603));
    set.add(new SimpleCnvRegion(807225, 807267));
    set.add(new SimpleCnvRegion(809048, 809078));
    set.add(new SimpleCnvRegion(809150, 809165));
    set.add(new SimpleCnvRegion(812588, 812631));
    final ComplexCnvRegion re = new ComplexCnvRegion(89, 1709001, set);
    re.globalIntegrity();

  }
}
