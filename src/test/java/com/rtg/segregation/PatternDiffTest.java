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
