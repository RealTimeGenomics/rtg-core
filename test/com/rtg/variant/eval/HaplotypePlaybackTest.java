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

package com.rtg.variant.eval;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class HaplotypePlaybackTest extends TestCase {

  public void testSimpleBase() {
    byte[] template = {1, 1, 1, 1, 1, 1, 2, 1, 1};
    HaplotypePlayback path = new HaplotypePlayback(template);
    //snp C at 2
    path.addVariant(new OrientedVariant(new MockVariant(1, 2, new byte[] {2}, null), true));
    //insert G at 4
    path.addVariant(new OrientedVariant(new MockVariant(4, 4, new byte[] {3}, null), true));
    //delete length 1 at 6
    path.addVariant(new OrientedVariant(new MockVariant(6, 7, new byte[] {}, null), true));

    byte[] expected = {2, 1, 1, 3, 1, 1, 2, 1, 1};
    check(expected, path);
  }

  public void testMoreComplex() {
    byte[] template = {1, 1, 1, 1, 1, 1, 2, 1, 1};
    HaplotypePlayback path = new HaplotypePlayback(template);
    //mnp A -> CTT:GGG at 2
    path.addVariant(new OrientedVariant(new MockVariant(1, 2, new byte[] {2, 4, 4}, new byte[] {3, 3, 3}), false));
    //insert A -> C:GA at 4
    path.addVariant(new OrientedVariant(new MockVariant(4, 5, new byte[] {2}, new byte[] {3, 1}), true));
    //delete A -> i:T
    path.addVariant(new OrientedVariant(new MockVariant(6, 7, new byte[] {}, new byte[] {4}), false));

    byte[] expected = {3, 3 , 3, 1, 1, 2, 1, 4, 2, 1, 1};
    check(expected, path);
  }
  public void testMoreComplexAlternate() {
                               // 1  2  3  4  5  6  7  8  9
    byte[] template = {1, 1, 1, 1, 1, 1, 2, 1, 1};
                               // 2441  1  31 4     3  1  1
    HaplotypePlayback path = new HaplotypePlayback(template);
    //mnp A -> CTT:GGG at 2
    path.addVariant(new OrientedVariant(new MockVariant(1, 2, new byte[] {2, 4, 4}, new byte[] {3, 3, 3}), true));
    //insert A -> C:GA at 4
    path.addVariant(new OrientedVariant(new MockVariant(4, 5, new byte[] {2}, new byte[] {3, 1}), false));

    //A -> G:T
    path.addVariant(new OrientedVariant(new MockVariant(5, 6, new byte[] {3}, new byte[] {4}), false));
    //delete A -> i:T
    path.addVariant(new OrientedVariant(new MockVariant(6, 7, new byte[] {}, new byte[] {4}), true));
    //C -> G:T
    path.addVariant(new OrientedVariant(new MockVariant(7, 8, new byte[] {3}, new byte[] {4}), true));

    byte[] expected = {2, 4 , 4, 1, 1, 3, 1, 4, 3, 1, 1};
    TestUtils.containsAll(path.toString()
        , "HaplotypePlayback: position=-1 inPosition=-1"
        , "current:1:2 CTT:GGG"
        , "future:");
    check(expected, path);
    TestUtils.containsAll(path.toString()
        , "HaplotypePlayback: position=8 inPosition=-1"
        , "current:" + null
        , "future:");
  }

  private void check(byte[] expected, HaplotypePlayback path) {
    int i = 0;
    while (path.hasNext()) {
      path.next();
      //System.err.println(path.nt());
      assertEquals(expected[i], path.nt());
      i++;
    }
    assertEquals(expected.length, i);
  }

  public void testCopy() {
    byte[] template = {1, 1, 1, 1, 1, 1, 2, 1, 1};
                               // 2441  1  31 1     3  1  1
    HaplotypePlayback path = new HaplotypePlayback(template);
    //mnp A -> CTT:GGG at 2
    path.addVariant(new OrientedVariant(new MockVariant(1, 2, new byte[] {2, 4, 4}, new byte[] {3, 3, 3}), true));
    assertEquals(path.toString(), path.copy().toString());

  }

  public void testPartialUpdated() {
                               // 1  2  3  4  5  6  7  8  9
    byte[] template = {1, 1, 1, 1, 1, 1, 2, 1, 1};
                               // 2441  1  31 1     3  1  1
    HaplotypePlayback path = new HaplotypePlayback(template);
    assertEquals(null, path.currentVariant());
    //mnp A -> CTT:GGG at 2
    path.addVariant(new OrientedVariant(new MockVariant(1, 2, new byte[] {2, 4, 4}, new byte[] {3, 3, 3}), true));
    //insert A -> C:GA at 4
    OrientedVariant first = new OrientedVariant(new MockVariant(4, 5, new byte[] {2}, new byte[] {3, 1}), false);
    path.addVariant(first);

    byte[] expected = {2, 4 , 4, 1, 1, 3, 1, 4, 3, 1, 1};
    int i = 0;
    while (i < 7) {
      path.next();
      //System.err.println(path.nt());
      assertEquals(expected[i], path.nt());
      i++;
    }
    assertEquals(first, path.currentVariant());
    //A -> G:T
    OrientedVariant next = new OrientedVariant(new MockVariant(5, 6, new byte[] {3}, new byte[] {4}), false);
    path.addVariant(next);
    //delete A -> i:T
    path.addVariant(new OrientedVariant(new MockVariant(6, 7, new byte[] {}, new byte[] {4}), true));
    //C -> G:T
    path.addVariant(new OrientedVariant(new MockVariant(7, 8, new byte[] {3}, new byte[] {4}), true));

    while (path.templatePosition() < 4) {
      path.next();
      //System.err.println(path.nt());
      assertEquals(expected[i], path.nt());
      i++;
    }
    assertEquals(next, path.currentVariant());

    while (path.hasNext()) {
      path.next();
      //System.err.println(path.nt());
      assertEquals(expected[i], path.nt());
      i++;
    }
    assertEquals(expected.length, i);
  }

  public void testComparable() {
                               // 1  2  3  4  5  6  7  8  9
    byte[] template = {1, 1, 1, 1, 1, 1, 2, 1, 1};
                               // 2441  1  31 1     3  1  1
    HaplotypePlayback path = new HaplotypePlayback(template);
    //mnp A -> CTT:GGG at 2
    path.addVariant(new OrientedVariant(new MockVariant(1, 2, new byte[] {2, 4, 4}, new byte[] {3, 3, 3}), true));
    //insert A -> C:GA at 4
    path.addVariant(new OrientedVariant(new MockVariant(4, 5, new byte[] {2}, new byte[] {3, 1}), false));

    assertEquals(0, path.compareTo(path));
    HaplotypePlayback copy = path.copy();
    assertTrue(path.equals(copy));
    assertTrue(path.hashCode() == copy.hashCode());
    copy.addVariant(new OrientedVariant(new MockVariant(5, 6, new byte[] {3}, new byte[] {4}), false));
    assertEquals(-1, path.compareTo(copy));
    assertEquals(1, copy.compareTo(path));
    assertFalse(path.equals(copy));
    assertFalse(path.hashCode() == copy.hashCode());
    assertFalse(path.equals(null));
    path.addVariant(new OrientedVariant(new MockVariant(5, 6, new byte[] {3}, new byte[] {4}), true));
    assertEquals(1, path.compareTo(copy));
    assertEquals(-1, copy.compareTo(path));

    HaplotypePlayback pathAdvanced = path.copy();
    pathAdvanced.next();
    assertEquals(-1, path.compareTo(pathAdvanced));
    assertEquals(1, pathAdvanced.compareTo(path));

    HaplotypePlayback path2 = new HaplotypePlayback(template);
    HaplotypePlayback copy2 = path2.copy();
    assertEquals(0, path2.compareTo(copy2));

    path2.addVariant(new OrientedVariant(new MockVariant(1, 2, new byte[] {2, 4, 4}, new byte[] {3, 3, 3}), true));
    assertEquals(1, path2.compareTo(copy2));
    assertEquals(-1, copy2.compareTo(path2));

    copy2.addVariant(new OrientedVariant(new MockVariant(2, 3, new byte[] {2, 4, 4}, new byte[] {3, 3, 3}), true));
    assertEquals(-1, path2.compareTo(copy2));
    assertEquals(1, copy2.compareTo(path2));

  }
  public void testComparable2() {
                               // 1  2  3  4  5  6  7  8  9
    byte[] template = {1, 1, 1, 1, 1, 1, 2, 1, 1};
                               // 2441  1  31 1     3  1  1
    HaplotypePlayback path = new HaplotypePlayback(template);
    path.addVariant(new OrientedVariant(new MockVariant(1, 2, new byte[] {2, 4, 4}, new byte[] {3, 3, 3}), true));
    HaplotypePlayback copy = path.copy();
    path.next();
    path.next();
    copy.next();
    assertEquals(1, path.compareTo(copy));
    assertEquals(-1, copy.compareTo(path));
  }

  public void testMoveForward() {
    byte[] template = {1, 1, 1, 1, 2, 3, 4, 1, 2, 3, 4, 1, 1, 2, 1, 1};
    HaplotypePlayback path = new HaplotypePlayback(template);
    //snp C at 0
    path.addVariant(new OrientedVariant(new MockVariant(1, 2, new byte[] {2}, null), true));
    //snp A at
    path.addVariant(new OrientedVariant(new MockVariant(14, 15, new byte[] {1}, null), true));

    byte[] expected = {2, 1, 1, 1, 2, 3, 4, 1, 2, 3, 4, 1, 1, 1, 1, 1};
    int i = -1;
    while (path.templatePosition() < 4) {
      path.next();
      i++;
      //System.err.println(path.nt());
      assertEquals("position: " + i, expected[i], path.nt());
      assertEquals(i, path.templatePosition());
    }
    i = 9;
    path.moveForward(i);
    while (path.hasNext()) {
      path.next();
      i++;
      assertEquals("position: " + i, expected[i], path.nt());
      assertEquals(i, path.templatePosition());
    }
    assertEquals(expected.length - 1, i);
  }
}
