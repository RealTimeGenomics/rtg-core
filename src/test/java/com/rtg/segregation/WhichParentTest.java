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
public class WhichParentTest extends TestCase {

  public void test() {
    final PatternDiff noDiffNoFlip = new PatternDiff(true, false, -1, false);
    final PatternDiff noDiffFlip = new PatternDiff(true, false, -1, true);
    final PatternDiff noExpl = new PatternDiff(false, true, -1, false);
    final PatternDiff validNoFlip = new PatternDiff(false, false, 42, false);
    final PatternDiff validFlip = new PatternDiff(false, false, 21, true);

    check(noDiffFlip, noDiffFlip, false, null);
    check(noDiffFlip, noExpl, false, null);
    check(noDiffFlip, validFlip, false, 21);
    check(noExpl, noExpl, false, null);
    check(noExpl, validFlip, false, null);
    check(validFlip, validFlip, false, null);

    check(noDiffNoFlip, noDiffFlip, false, null);
    check(noDiffNoFlip, noExpl, false, null);
    check(noDiffNoFlip, validFlip, false, 21);
    check(noDiffNoFlip, validNoFlip, false, 42);
    check(noExpl, noExpl, false, null);
    check(noExpl, validNoFlip, false, null);
    check(validFlip, validNoFlip, false, null);

  }

  private void check(final PatternDiff a, final PatternDiff b, boolean expFather, Integer expChild) {
    checkx(a, b, expFather, expChild);
    checkx(b, a, !expFather, expChild);
  }

  private void checkx(final PatternDiff a, final PatternDiff b, boolean expFather, Integer expChild) {
    final WhichParent t0 = new WhichParent(a, b);
    t0.integrity();
    assertTrue(a == t0.father());
    assertTrue(b == t0.mother());
    if (expChild == null) {
      assertFalse(t0.isValid());
      return;
    }
    assertTrue(t0.isValid());
    assertEquals(expFather, t0.isFather());
    assertEquals(expChild.intValue(), t0.child());
  }
}
