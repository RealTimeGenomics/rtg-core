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

package com.rtg.segregation;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class PatternArrayTest extends TestCase {

  public void test() {
    final PatternArray pa = new PatternArray(new Pattern(8), new Pattern(2));
    assertEquals(2, pa.length());
    assertEquals(8, pa.index(0).intPattern());
    assertEquals(2, pa.index(1).intPattern());
    assertTrue(pa.isUniqueAndValid(false));
    assertEquals(" fa: 10 mo: 11", pa.toString());
  }

  //Test String based constructor
  public void testStr() {
    final PatternArray pa = new PatternArray(new Pattern(15), new Pattern(2));
    assertEquals(2, pa.length());
    assertEquals(15, pa.index(0).intPattern());
    assertEquals(2, pa.index(1).intPattern());
  }

  public void testFaString() {
    final PatternArray pa = new PatternArray(new Pattern(15), new Pattern(2));
    assertEquals("?0", pa.faString());
    assertFalse(pa.isUniqueAndValid(false));
    assertEquals(1, pa.faCount());
    assertEquals(null, pa.phaseCount(true));
  }

  public void testFaCount() {
    final PatternArray pa = new PatternArray(new Pattern(2), new Pattern(2));
    assertEquals("00", pa.faString());
    assertEquals(2, pa.faCount());
  }

  public void testMoString() {
    final PatternArray pa = new PatternArray(new Pattern(15), new Pattern(2));
    assertEquals("?1", pa.moString());
    assertFalse(pa.isUniqueAndValid(false));
    assertEquals(0, pa.moCount());
    assertEquals(null, pa.phaseCount(false));
  }

  public void testMoCount() {
    final PatternArray pa = new PatternArray(new Pattern(1), new Pattern(1));
    assertEquals("00", pa.moString());
    assertEquals(2, pa.moCount());
  }

  public void testPhaseCount() {
    final PatternArray pa = new PatternArray(new Pattern(2), new Pattern(2));
    assertEquals(Integer.valueOf(2), pa.phaseCount(true));
    assertEquals(Integer.valueOf(0), pa.phaseCount(false));
  }

  public void testCompatible() {
    checkCompatible(new int[] {0b1111}, new int[] {0b1111}, true);
    checkCompatible(new int[] {0b1000}, new int[] {0b1000}, true);
    checkCompatible(new int[] {0b1000, 0b1000, 0b1000, 0b1000}, new int[] {0b0100, 0b0100, 0b0100, 0b0100}, true);
    checkCompatible(new int[] {0b1000, 0b1000, 0b1000, 0b1000}, new int[] {0b0100, 0b0010, 0b0100, 0b0100}, false);
    checkCompatible(new int[] {0b1000, 0b1000, 0b1000, 0b1000}, new int[] {0b0100, 0b0100, 0b0100, 0b0010}, false);
    checkCompatible(new int[] {0b1000, 0b1000, 0b1000, 0b1000}, new int[] {0b0010, 0b0100, 0b0100, 0b0100}, false);
    checkCompatible(new int[] {0b0100, 0b0100, 0b0100, 0b0100}, new int[] {0b0100, 0b0100, 0b0100, 0b0100}, true);
    checkCompatible(new int[] {0b0100, 0b0100, 0b0100, 0b0100}, new int[] {0b1000, 0b0100, 0b0100, 0b0100}, false);
  }

  private void checkCompatible(int[] a, int[] b, boolean exp) {
    final PatternArray pa = pat(a);
    final PatternArray pb = pat(b);
    assertEquals(exp, pa.compatible(pb));
    assertEquals(exp, pb.compatible(pa));
  }

  PatternArray pat(final int[] q) {
    final Pattern[] p = new Pattern[q.length];
    for (int i = 0; i < q.length; ++i) {
      p[i] = new Pattern(q[i]);
    }
    return new PatternArray(p);
  }

  public void testFlipIntersect() {
    checkFlipIntersect(new int[] {15, 15, 15, 15}, new int[] {1, 2, 4, 8}, 0, new int[] {1, 2, 4, 8});
    checkFlipIntersect(new int[] {15, 15, 15, 15}, new int[] {1, 2, 4, 8}, 1, new int[] {4, 8, 1, 2});
    checkFlipIntersect(new int[] {15, 15, 15, 15}, new int[] {1, 2, 4, 8}, 2, new int[] {2, 1, 8, 4});
    checkFlipIntersect(new int[] {15, 15, 15, 15}, new int[] {1, 2, 4, 8}, 3, new int[] {8, 4, 2, 1});
    checkFlipIntersect(new int[] {0, 0}, new int[] {0, 15}, 0, null);
  }

  void checkFlipIntersect(final int[] a, final int[] b, final int f, final int[] exp) {
    final PatternArray ap = int2PA(a);
    final PatternArray bp = int2PA(b);
    final PatternArray ex = int2PA(exp);
    final PatternArray ac = ap.flipIntersect(bp, f);
    if (ac != null) {
      if (exp == null) {
        fail();
      } else {
        assertTrue(ex.strictEquals(ac));
      }
    }
  }

  private PatternArray int2PA(final int[] q) {
    if (q == null) {
      return null;
    }
    return new PatternArray(q);
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
    checkFlipIntersect(new int[] {a}, new int[] {b}, f, new int[] {exp});
  }

  public void testEqualsHash() {
    checkFlipEquals(1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
    TestUtils.equalsHashTest(new PatternArray[][] {new PatternArray[] {new PatternArray(1, 2, 3)}, new PatternArray[] {new PatternArray(7, 9, 11)}});
  }

  private void checkFlipEquals(int... q) {
    final PatternArray pa = new PatternArray(q);
    assertFalse(pa.strictEquals(null));
    assertTrue(pa.strictEquals(pa));
    for (int f = 0; f < Pattern.NUMBER_FLIPS; ++f) {
      final PatternArray pf = pa.flip(f);
      assertEquals(pa, pf);
      assertEquals(pf, pa);
      assertEquals(pa.hashCode(), pf.hashCode());
      if (f != 0) {
        assertFalse(pa.strictEquals(pf));
        assertFalse(pf.strictEquals(pa));
      }
    }
  }

  public void testHash() {
    assertEquals(685811806, new PatternArray(1, 2, 3).hashCode());
  }

  public void testStringToPattern() {
    checkStringToPattern(1, 2, 4, 8);
  }

  private void checkStringToPattern(int... q) {
    final PatternArray pa = new PatternArray(q);
    final StringBuilder fa = new StringBuilder();
    final StringBuilder mo = new StringBuilder();
    for (int i = 0; i < q.length; ++i) {
      fa.append(pa.index(i).faString());
      mo.append(pa.index(i).moString());
    }
    final PatternArray pb = new PatternArray(fa.toString(), mo.toString());
    assertTrue(pa.strictEquals(pb));
  }

  public void testCrossover() {
    checkXO("00100001110", "10000110100", "00100001100", "11100110100", null, null);
    checkXO("00100001110", "10000110100", "00100001110", "10000110100", null, null); //equal
    checkXO("00100001110", "10000110100", "00110001110", "10000110100", "fa 3 7 6", "fa 3 4 5");
    checkXO("00100001110", "10000110100", "00100001110", "11000110100", "mo 1 7 6", "mo 1 4 5");
    checkXO("00100001110", "10000110100", "10100001110", "10000110100", "fa 0 7 6", "fa 0 4 5");
    checkXO("00100001110", "10000110100", "00100001110", "00000110100", "mo 0 7 8", "mo 0 4 3");
  }

  private void checkXO(final String fa, final String ma, final String fb, final String mb, final String exp0, final String exp1) {
    final String faf = flip(fa);
    final String maf = flip(ma);
    final String fbf = flip(fb);
    final String mbf = flip(mb);
    checkXO0(fa, ma, fb, mb, exp0, exp1);
    checkXO0(fa, ma, fb, mbf, exp0, exp1);
    checkXO0(fa, ma, fbf, mb, exp0, exp1);
    checkXO0(fa, ma, fbf, mbf, exp0, exp1);
    checkXO0(fa, maf, fb, mb, exp0, exp1);
    checkXO0(fa, maf, fb, mbf, exp0, exp1);
    checkXO0(fa, maf, fbf, mb, exp0, exp1);
    checkXO0(fa, maf, fbf, mbf, exp0, exp1);
    checkXO0(faf, ma, fb, mbf, exp0, exp1);
    checkXO0(faf, ma, fbf, mb, exp0, exp1);
    checkXO0(faf, ma, fbf, mbf, exp0, exp1);
    checkXO0(faf, maf, fb, mb, exp0, exp1);
    checkXO0(faf, maf, fb, mbf, exp0, exp1);
    checkXO0(faf, maf, fbf, mb, exp0, exp1);
    checkXO0(faf, maf, fbf, mbf, exp0, exp1);
  }

  String flip(final String a) {
    return a.replace('0', 'X').replace('1', '0').replace('X', '1');
  }

  private void checkXO0(final String fa, final String ma, final String fb, final String mb, final String exp0, final String exp1) {
    final PatternArray a = new PatternArray(fa, ma);
    final PatternArray b = new PatternArray(fb, mb);
    final CrossOver xo = PatternArray.crossover(a, b, false);
    if (exp0 == null) {
      assertNull(xo);
    } else {
      final String xos = xo.toString();
      assertTrue(xos, exp0.equals(xos) || exp1.equals(xos));
    }
  }

  public void testCrossOverInvalid() {
    checkCrossOverInvalid(2, 4, true);
    checkCrossOverInvalid(2, 15, true);
  }

  private void checkCrossOverInvalid(final int a, final int b, final boolean isNull) {
    final PatternArray pa = new PatternArray(new Pattern(a));
    final PatternArray pb = new PatternArray(new Pattern(b));
    final CrossOver xoa = PatternArray.crossover(pa, pb, false);
    assertEquals(isNull, xoa == null);
    final CrossOver xob = PatternArray.crossover(pb, pa, false);
    assertEquals(isNull, xob == null);
  }

  //seen in a real case
  public void testCrossOverBug() {
    final PatternArray a = new PatternArray("00100001110", "10000110100");
    final PatternArray b = new PatternArray("00100001100", "11100110100");
    final CrossOver xo = PatternArray.crossover(a, b, false);
    assertNull(xo);
  }

  public void testParital() {
    final PatternArray a = new PatternArray("000???0?0??", "00110001101");
    assertEquals("000???0?0??", a.faString());
    assertEquals("00110001101", a.moString());
    assertTrue(a.isUniqueAndValid(true));
  }
}
