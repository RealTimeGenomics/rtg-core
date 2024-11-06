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

package com.rtg.variant.sv;

import com.rtg.util.integrity.Exam;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractSamCountsTest extends TestCase {

  protected abstract SamCounts getCounts(int length);

  public void test0() {
    final SamCounts sa = getCounts(0);
    Exam.integrity(sa);
    assertEquals(0.0, sa.count(0, 0));

    try {
      sa.increment(0);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("index=0 length=0", e.getMessage());
    }
  }

  public void test() {
    final SamCounts sa = getCounts(10);
    Exam.integrity(sa);
    assertEquals(0.0, sa.count(5, 0));
    assertEquals(0.0, sa.count(0, -1));
    assertEquals(0.0, sa.count(0, 0));
    assertEquals(0.0, sa.count(0, 0));
    assertEquals(0.0, sa.count(0, 0));
    assertEquals(0.0, sa.count(0, 9));
    assertEquals(0.0, sa.count(0, 10));

    try {
      sa.increment(-1);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("index=-1 length=10", e.getMessage());
    }

    try {
      sa.increment(10);
      fail();
    } catch (final IndexOutOfBoundsException e) {
      assertEquals("index=10 length=10", e.getMessage());
    }

    check10(sa);
  }

  public void testFlushTo() {
    checkFlush(10, 5, 5);
    checkFlush(10, 7, 5);
  }

  private void checkFlush(int max, int buffer, int flushpoint) {
    final SamCounts sa = getCounts(0);
    sa.reset(max, buffer);
    for (int i = 0; i < flushpoint; ++i) {
      sa.increment(i);
      assertEquals(1.0, sa.count(i, 0));
    }
    assertEquals(1.0, sa.count(flushpoint - 1, 0)); // Already flushed
    assertEquals(0.0, sa.count(flushpoint, 0));
    sa.flushTo(flushpoint);
    assertEquals(0.0, sa.count(flushpoint - 1, 0)); // Already flushed
    assertEquals(0.0, sa.count(flushpoint, 0)); // Should be OK, nothing inserted yet

    for (int i = flushpoint; i < max; ++i) {
      sa.increment(i);
      assertEquals(1.0, sa.count(i, 0));
    }
    assertEquals(0.0, sa.count(flushpoint - 1, 0)); // Already flushed
  }


  public void testReset1() {
    final SamCounts sa = getCounts(2);
    Exam.integrity(sa);

    sa.reset(10, 10);
    check10(sa);
  }

  public void testReset2() {
    final SamCounts sa = getCounts(20);
    Exam.integrity(sa);

    sa.reset(10, 10);
    check10(sa);
  }

  protected void check10Inc(final SamCounts sa) {
    sa.increment(0);
    sa.increment(4);
    sa.increment(4);
  }
  protected void check10(final SamCounts sa) {
    check10Inc(sa);

    assertEquals(1.0, sa.count(0, 0));
    assertEquals(2.0, sa.count(4, 0));
    assertEquals(1.0, sa.count(4, -4, 0));
    assertEquals(3.0, sa.count(4, -4, 1));
    assertEquals(2.0, sa.count(3, 1));
    assertEquals(2.0, sa.count(5, -1));
    assertEquals(0.0, sa.count(5, 0));
    assertEquals(0.0, sa.count(4, 1));
    assertEquals(0.0, sa.count(6, 0));
    assertEquals(0.0, sa.count(9, 0));
    assertEquals(0.0, sa.count(15, 0));

    assertEquals(0.0, sa.sumLn(4, -4, 0), 0.00001);
    assertEquals(1.38629436, sa.sumLn(5, -1, 1), 0.00001);
  }

}
