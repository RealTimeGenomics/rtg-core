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
package com.rtg.metagenomics;

import java.util.ArrayList;
import java.util.Arrays;

import com.rtg.metagenomics.matrix.MatrixSymmetric;
import com.rtg.metagenomics.matrix.Vector;
import com.rtg.util.SortedMultiSet;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class FragTest extends TestCase {

  public void test() {
    final Frag fr = frag(1, 1, 42);
    assertEquals("  1:2  42:1 {1}", fr.toString());
    fr.checkRange(43);
    fr.globalIntegrity();
    try {
      fr.checkRange(42);
      fail();
    } catch (final Exception e) {
      //expected
      assertEquals(null, e.getMessage());
    }
    assertEquals(1, fr.multiplicity());
    fr.setMultiplicity(42);
    assertEquals(42, fr.multiplicity());
    @SuppressWarnings("unchecked")
    final SortedMultiSet<Integer>[] counts = (SortedMultiSet<Integer>[]) new SortedMultiSet<?>[43];
    for (int i = 0; i < counts.length; ++i) {
      counts[i] = new SortedMultiSet<>();
    }
    fr.count(counts);
    assertEquals("[]", counts[0].toString());
    assertEquals("[2->42]", counts[1].toString());
    assertEquals("[2->42]", counts[42].toString());
    final Frag fr0 = frag(1, 1, 42);
    final Frag fr1 = frag(1, 2, 42);
    final Frag fr2 = frag(42);
    TestUtils.equalsTest(new Frag[][] {{fr, fr0}, {fr1}, {fr2}});
  }

  static Frag frag(final int...f) {
    final ArrayList<Integer> al = new ArrayList<>();
    for (final int fi : f) {
      al.add(fi);
    }
    return new Frag(al);
  }

  public void testStatsSimple() {
    final Frag fr = frag(1, 1, 42);
    final int[] a = new int[43];
    Arrays.fill(a, -1);
    fr.stats(a);
    assertEquals(-1, a[0]);
    assertEquals(1, a[1]);
    assertEquals(1, a[42]);
    for (int i = 0; i < a.length; ++i) {
      assertTrue(a[i] <= i);
    }
    //System.err.println(Arrays.toString(a));
  }


  public void testStatsComplicated() {
    final int[] a = new int[7];
    Arrays.fill(a, -1);
    frag(0, 1, 2).stats(a);
    checkA(a, 0, 0, 0, -1, -1, -1, -1);
    frag(3, 4).stats(a);
    checkA(a, 0, 0, 0, 3, 3, -1, -1);
    frag(5, 6).stats(a);
    checkA(a, 0, 0, 0, 3, 3, 5, 5);
    frag(3, 6).stats(a);
    checkA(a, 0, 0, 0, 3, 3, 3, 3);
    frag(1, 5).stats(a);
    checkA(a, 0, 0, 0, 0, 3, 0, 3);
    frag(4, 5).stats(a);
    checkA(a, 0, 0, 0, 0, 0, 0, 3);
  }


  public void testStatsSingle() {
    final int[] a = new int[2];
    Arrays.fill(a, -1);
    frag(0).stats(a);
    checkA(a, 0, -1);
    frag(1).stats(a);
    checkA(a, 0, 1);
    frag(1, 0).stats(a);
    checkA(a, 0, 0);
  }

  private void checkA(final int[] a, final int... b) {
    assertEquals(b.length, a.length);
    for (int i = 0; i < a.length; ++i) {
      assertEquals("" + i, a[i], b[i]);
    }
    for (int i = 1; i < a.length; ++i) {
      assertTrue("" + i, a[i] <= i);
    }
  }

  public void testIncrement() {
    final ArrayList<Integer> al = new ArrayList<>();
    al.add(1);
    al.add(1);
    al.add(42);
    final Frag fr = new Frag(al);
    final Vector j = new Vector(44);
    final Vector r = new Vector(44);
    r.set(1, 0.5);
    r.set(42, 3.0);
    assertEquals(-Math.log(4.0), fr.increment(r, j));
    assertEquals(-Math.log(4.0), fr.l(r));
    assertEquals(0.0, j.get(0));
    assertEquals(-0.5, j.get(1));
    assertEquals(-0.25, j.get(42));
  }

  public void testSum() {
    final ArrayList<Integer> al = new ArrayList<>();
    al.add(1);
    al.add(1);
    al.add(42);
    final Frag fr = new Frag(al);
    final Vector j = new Vector(44);
    j.set(1, 0.75);
    j.set(42, 0.1);
    assertEquals(1.6, fr.sum(j));
  }

  public void testEquals() {
    final ArrayList<Integer> al = new ArrayList<>();
    al.add(1);
    al.add(1);
    al.add(42);
    final Frag fr1 = new Frag(al);
    final Frag fr2 = new Frag(al);
    assertEquals(fr1, fr2);
    assertEquals(fr1.hashCode(), fr2.hashCode());
    al.add(43);
    final Frag fr3 = new Frag(al);
    assertFalse(fr1.equals(fr3));
    assertTrue(fr1.hashCode() != fr3.hashCode());
  }

  public void testIncrementWithMultiplicity() {
    final ArrayList<Integer> al = new ArrayList<>();
    al.add(1);
    al.add(1);
    al.add(42);
    final Frag fr1 = new Frag(al);
    final Frag fr2 = new Frag(al);
    final Frag fr3 = new Frag(al);
    fr3.setMultiplicity(2);
    final Vector j1 = new Vector(44);
    final Vector r1 = new Vector(44);
    r1.set(1, 0.5);
    r1.set(42, 3.0);
    assertEquals(-Math.log(4.0), fr1.increment(r1, j1));
    assertEquals(-Math.log(4.0), fr2.increment(r1, j1));
    assertEquals(-Math.log(4.0), fr2.l(r1));
    assertEquals(0.0, j1.get(0));
    assertEquals(-1.0, j1.get(1));
    assertEquals(-0.5, j1.get(42));
    final Vector j2 = new Vector(44);
    final Vector r2 = new Vector(44);
    r2.set(1, 0.5);
    r2.set(42, 3.0);
    assertEquals(-2 * Math.log(4.0), fr3.increment(r2, j2));
    assertEquals(-2 * Math.log(4.0), fr3.l(r2));
    assertEquals(j1.get(0), j2.get(0));
    assertEquals(j1.get(1), j2.get(1));
    assertEquals(j1.get(42), j2.get(42));
    assertEquals(r1.get(1), r2.get(1));
    assertEquals(r1.get(42), r2.get(42));

    fr1.init(j1);
    fr2.init(j1);
    fr3.init(j2);
    assertEquals(j1.get(0), j2.get(0));
    assertEquals(j1.get(1), j2.get(1));
    assertEquals(j1.get(42), j2.get(42));

    final MatrixSymmetric h1 = new MatrixSymmetric(44);
    final MatrixSymmetric h2 = new MatrixSymmetric(44);
    // Just a nonzero start position, symmetric
    for (int k = 0; k < 44; ++k) {
      for (int j = 0; j <= k; ++j) {
        h1.set(k, j, (0.5 + k + j) / (0.5 + k * j));
        h2.set(k, j, (0.5 + k + j) / (0.5 + k * j));
      }
    }
    assertEquals(-Math.log(4.0), fr1.increment(r1, j1, h1));
    assertEquals(-Math.log(4.0), fr2.increment(r1, j1, h1));
    assertEquals(-2 * Math.log(4.0), fr3.increment(r1, j1, h2));
    for (int k = 0; k < 44; ++k) {
      for (int j = 0; j < 44; ++j) {
        assertEquals(h1.get(k, j), h2.get(k, j), 1e-6);
      }
    }
  }
  public void testSubBlock() {
    final Frag f = frag(0, 2, 4);
    assertEquals(0, f.subBlock(0, 41, 0, 42, 0));
    assertEquals(1, f.subBlock(1, 41, 1, 42, 1));
  }
  public void testSubBlockNotFirst() {
    final Frag f = frag(1, 2, 4);
    assertEquals(41, f.subBlock(0, 41, 41, 42, 41));
  }
  public void testSubBlockLengthOne() {
    final Frag f = frag(0);
    assertEquals(0, f.subBlock(0, 41, 0, 42, 0));
    assertEquals(1, f.subBlock(1, 41, 1, 42, 1));
  }

  public void testSubFrag() {
    final Frag f = frag(0, 0, 2, 4);
    f.setMultiplicity(42);
    final Frag g = frag(0, 0, 1, 2);
    g.setMultiplicity(42);
    final Frag actual = f.subFrag(0, 101, 1, 101, 2);
    assertTrue(g.identical(actual));
  }
}
