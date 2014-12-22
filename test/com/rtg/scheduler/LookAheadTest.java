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

package com.rtg.scheduler;

import junit.framework.TestCase;

/**
 */
public class LookAheadTest extends TestCase {

  public void test() {
    final LookAhead la = new LookAhead(3, 0);
    assertEquals(3, la.lookAhead());
    checkEmpty(la);
    assertTrue(la.ok(0, 0));
    assertTrue(la.ok(3, 0));

    la.increment(0);
    check(la, 0, 1);

    la.increment(1);
    check(la, 0, 2);

    la.increment(2);
    check(la, 0, 3);

    la.decrement(0);
    check(la, 1, 2);

    la.decrement(1);
    check(la, 2, 1);

    la.decrement(2);
    checkEmpty(la);

    la.increment(3);
    //System.err.println(la.toString());
    check(la, 3, 1);
    try {
      la.increment(2);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Time invalid: current=3 t=2 lookAhead=3 mDelta=0", e.getMessage());
    }
    try {
      la.increment(7);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Time invalid: current=3 t=7 lookAhead=3 mDelta=0", e.getMessage());
    }
    la.increment(4);
    check(la, 3, 2);
    la.increment(5);
    check(la, 3, 3);
    la.decrement(3);
    check(la, 4, 2);

    la.increment(6);
    check(la, 4, 3);

    la.decrement(5);
    check(la, 4, 2);
    la.decrement(4);
    check(la, 6, 1);
  }

  public void testDelta() {
    final int delta = 2;
    final LookAhead la = new LookAhead(3, delta);
    assertEquals(3, la.lookAhead());
    assertEquals(delta, la.delta());
    checkEmpty(la);

    la.increment(3);
    //System.err.println(la.toString());
    check(la, 3, 1);
    final int badTime = 3 + la.lookAhead() + la.delta() + 1;
    try {
      la.increment(badTime);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Time invalid: current=3 t=" + badTime + " lookAhead=3 mDelta=2", e.getMessage());
    }
    la.increment(8);
    check(la, 3, 2);
  }

  //set up a long gap so check the loop for moving forward
  public void testLong() {
    final LookAhead la = new LookAhead(10, 0);
    assertEquals(10, la.lookAhead());

    for (int i = 0; i < 10; i++) {
      la.increment(i);
      la.globalIntegrity();
    }
    check(la, 0, 10);

    for (int i = 7, c = 9; i > 0; i--, c--) {
      la.decrement(i);
      check(la, 0, c);
    }
    la.decrement(0);
    check(la, 8, 2);
  }

  //set up a long gap so check the loop for moving forward
  public void testLong1() {
    final LookAhead la = new LookAhead(10, 0);
    assertEquals(10, la.lookAhead());

    for (int i = 0; i < 10; i++) {
      la.increment(i);
      la.globalIntegrity();
    }
    check(la, 0, 10);

    for (int i = 8, c = 9; i > 0; i--, c--) {
      la.decrement(i);
      check(la, 0, c);
    }
    la.decrement(0);
    check(la, 9, 1);
  }

  //set up a long gap so check the loop for moving forward
  public void testLong10() {
    final LookAhead la = new LookAhead(10, 0);
    assertEquals(10, la.lookAhead());

    for (int i = 0; i <= 10; i++) {
      la.increment(i);
      la.globalIntegrity();
    }
    check(la, 0, 11);

    for (int i = 9, c = 10; i > 0; i--, c--) {
      la.decrement(i);
      check(la, 0, c);
    }
    la.decrement(0);
    check(la, 10, 1);
  }

  //0 look ahead
  public void test0() {
    final LookAhead la = new LookAhead(0, 0);
    assertEquals(0, la.lookAhead());
    checkEmpty(la);
    assertTrue(la.ok(0, 0));
    assertTrue(la.ok(3, 0));

    la.increment(0);
    check(la, 0, 1);

    la.decrement(0);
    checkEmpty(la);

    la.increment(1);
    check(la, 1, 1);

    la.decrement(1);
    checkEmpty(la);
  }

  private void checkEmpty(final LookAhead la) {
    la.globalIntegrity();
    assertEquals(0, la.total());
  }

  private void check(final LookAhead la, final int t, final int cnt) {
    la.globalIntegrity();
    assertEquals(cnt, la.total());
    if (t >= 1) {
      assertFalse(la.ok(t - 1, la.delta()));
    }
    assertTrue(la.ok(t, la.delta()));
    assertTrue(la.ok(t + la.lookAhead(), la.delta()));
    assertTrue(la.ok(t + la.lookAhead() + la.delta(), la.delta()));
    assertFalse(la.ok(t + la.lookAhead() + la.delta() + 1, la.delta()));
  }

}
