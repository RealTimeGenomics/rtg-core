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
