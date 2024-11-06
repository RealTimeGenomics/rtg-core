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

    for (int i = 0; i < 10; ++i) {
      la.increment(i);
      la.globalIntegrity();
    }
    check(la, 0, 10);

    for (int i = 7, c = 9; i > 0; --i, --c) {
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

    for (int i = 0; i < 10; ++i) {
      la.increment(i);
      la.globalIntegrity();
    }
    check(la, 0, 10);

    for (int i = 8, c = 9; i > 0; --i, --c) {
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

    for (int i = 0; i <= 10; ++i) {
      la.increment(i);
      la.globalIntegrity();
    }
    check(la, 0, 11);

    for (int i = 9, c = 10; i > 0; --i, --c) {
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
