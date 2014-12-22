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

package com.rtg.variant.realign.treealign;

import com.rtg.util.TestUtils;


/**
 */
public class SubTrieTest extends TrieTest {

  public void testSubHashCode() {
    final SubTrie st0 = new SubTrie(2, 1);
    assertEquals(System.identityHashCode(st0), st0.hashCode());
  }

  public void testSubEquality1() {
    final SubTrie st0 = new SubTrie(2, 1);
    final SubTrie st1 = new SubTrie(3, 1);
    final SubTrie st2 = new SubTrie(2, 2);
    TestUtils.equalsTest(new SubTrie[] {st0, st1, st2});
  }

  //differences at the level of children pointers
  public void testSubEqualityCh0() {
    final SubTrie st0 = new SubTrie(2, 1);
    final SubTrie st1 = new SubTrie(2, 2);
    final SubTrie st2 = new SubTrie(2, 2);
    st1.mChildren[1] = st0;
    assertFalse(st1.equals(st2));
    assertFalse(st2.equals(st1));
  }

  //the children are different
  public void testSubEqualityCh1a() {
    final SubTrie st0a = new SubTrie(2, 1);
    final SubTrie st0b = new SubTrie(1, 1);
    st0a.mChildren[3] = st0b;

    final SubTrie st1a = new SubTrie(2, 1);
    final SubTrie st1b = new SubTrie(1, 0);
    st1a.mChildren[3] = st1b;

    assertFalse(st1a.equals(st0a));
    assertFalse(st0a.equals(st1a));
  }

  //the children are different
  public void testSubEqualityCh1b() {
    final SubTrie st0a = new SubTrie(2, 1);
    final SubTrie st0b = new SubTrie(1, 1);
    st0a.mChildren[0] = st0b;

    final SubTrie st1a = new SubTrie(2, 1);
    final SubTrie st1b = new SubTrie(1, 0);
    st1a.mChildren[0] = st1b;

    assertFalse(st1a.equals(st0a));
    assertFalse(st0a.equals(st1a));
  }

  //the children are the same
  public void testSubEqualityCh2() {
    final SubTrie st0a = new SubTrie(2, 1);
    final SubTrie st0b = new SubTrie(1, 1);
    st0a.mChildren[3] = st0b;

    final SubTrie st1a = new SubTrie(2, 1);
    final SubTrie st1b = new SubTrie(1, 1);
    st1a.mChildren[3] = st1b;

    assertTrue(st1a.equals(st0a));
    assertTrue(st0a.equals(st1a));
  }

  public void testSubEqualityNull() {
    final SubTrie st0 = new SubTrie(2, 1);
    assertFalse(st0.equals(null));
    assertTrue(st0.equals(st0));
  }
}
