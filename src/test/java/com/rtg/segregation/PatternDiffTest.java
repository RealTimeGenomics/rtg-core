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

import junit.framework.TestCase;

/**
 */
public class PatternDiffTest extends TestCase {

  public void test() {
    final PatternDiff pd = new PatternDiff(false, false, 42, true);
    pd.integrity();
    assertFalse(pd.noDifference());
    assertFalse(pd.noExplantion());
    assertTrue(pd.isValid());
    assertTrue(pd.flipped());
    assertEquals(42, pd.child());
  }

  public void testNoDifference() {
    final PatternDiff pd = new PatternDiff(true, false, -1, true);
    pd.integrity();
    assertTrue(pd.noDifference());
    assertFalse(pd.noExplantion());
    assertFalse(pd.isValid());
    assertTrue(pd.flipped());
    try {
      pd.child();
      fail();
    } catch (final RuntimeException e) {
      assertEquals("No child defined", e.getMessage());
    }
  }

  public void testNoExplantion() {
    final PatternDiff pd = new PatternDiff(false, true, -1, false);
    pd.integrity();
    assertFalse(pd.noDifference());
    assertTrue(pd.noExplantion());
    assertFalse(pd.isValid());
    try {
      pd.child();
      fail();
    } catch (final RuntimeException e) {
      assertEquals("No child defined", e.getMessage());
    }
    try {
      pd.flipped();
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Flip not defined", e.getMessage());
    }
  }
}
