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
    assertFalse(re.contains(Integer.MIN_VALUE));
    assertFalse(re.contains(-1));
    assertTrue(re.contains(0));
    assertFalse(re.contains(1));
    assertFalse(re.contains(2));
    assertFalse(re.contains(3));
    assertFalse(re.contains(4));
    assertFalse(re.contains(5));
    assertFalse(re.contains(6));
    assertTrue(re.contains(7));
    assertTrue(re.contains(8));
    assertFalse(re.contains(9));
    assertFalse(re.contains(10));
    assertFalse(re.contains(11));
    assertFalse(re.contains(12));
    assertFalse(re.contains(13));
    assertFalse(re.contains(Integer.MAX_VALUE));

    assertEquals(0, re.getStart());
    assertEquals(12, re.getEnd());

    for (int i = 0; i <= 12; ++i) {
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
    assertFalse(re.contains(Integer.MIN_VALUE));
    assertFalse(re.contains(-1));
    assertTrue(re.contains(0));
    assertFalse(re.contains(1));
    assertFalse(re.contains(2));
    assertTrue(re.contains(3));
    assertTrue(re.contains(4));
    assertFalse(re.contains(5));
    assertFalse(re.contains(6));
    assertTrue(re.contains(7));
    assertTrue(re.contains(8));
    assertFalse(re.contains(9));
    assertFalse(re.contains(10));
    assertFalse(re.contains(11));
    assertFalse(re.contains(12));
    assertFalse(re.contains(13));
    assertFalse(re.contains(Integer.MAX_VALUE));

    assertEquals(0, re.getStart());
    assertEquals(12, re.getEnd());

    for (int i = 0; i <= 12; ++i) {
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
    assertFalse(re.contains(Integer.MIN_VALUE));
    assertFalse(re.contains(-5));
    assertFalse(re.contains(-4));
    assertFalse(re.contains(-3));
    assertFalse(re.contains(-2));
    assertFalse(re.contains(-1));
    assertFalse(re.contains(0));
    assertTrue(re.contains(1));
    assertFalse(re.contains(2));
    assertTrue(re.contains(3));
    assertTrue(re.contains(4));
    assertTrue(re.contains(5));
    assertTrue(re.contains(6));
    assertTrue(re.contains(7));
    assertTrue(re.contains(8));
    assertFalse(re.contains(9));
    assertFalse(re.contains(10));
    assertTrue(re.contains(11));
    assertTrue(re.contains(12));
    assertFalse(re.contains(13));
    assertFalse(re.contains(14));
    assertFalse(re.contains(15));
    assertFalse(re.contains(16));
    assertFalse(re.contains(17));
    assertFalse(re.contains(18));
    assertFalse(re.contains(Integer.MAX_VALUE));

    assertEquals(1, re.getStart());
    assertEquals(12, re.getEnd());

    for (int i = 2; i <= 12; ++i) {
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
    for (int i = 0; i < 2640; ++i) {
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
