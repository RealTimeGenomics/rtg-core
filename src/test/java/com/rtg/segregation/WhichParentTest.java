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
