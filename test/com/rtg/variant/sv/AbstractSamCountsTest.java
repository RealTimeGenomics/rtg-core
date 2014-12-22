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
    for (int i = 0; i < flushpoint; i++) {
      sa.increment(i);
      assertEquals(1.0, sa.count(i, 0));
    }
    assertEquals(1.0, sa.count(flushpoint - 1, 0)); // Already flushed
    assertEquals(0.0, sa.count(flushpoint, 0));
    sa.flushTo(flushpoint);
    assertEquals(0.0, sa.count(flushpoint - 1, 0)); // Already flushed
    assertEquals(0.0, sa.count(flushpoint, 0)); // Should be OK, nothing inserted yet

    for (int i = flushpoint; i < max; i++) {
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
