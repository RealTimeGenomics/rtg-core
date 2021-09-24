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
public class ResultTest extends TestCase {

  public void test1null() {
    final Result r = new Result((Object) null);
    assertEquals(null, r.result(0));
    assertEquals(1, r.length());
  }

  public void test1() {
    final Result r = new Result(42);
    assertEquals(42, r.result(0));
    assertEquals(1, r.length());
    assertEquals("[42]", r.toString());
  }

  public void test0() {
    final Result r = new Result();
    assertEquals(0, r.length());
    try {
      r.result(0);
      fail();
    } catch (final RuntimeException e) {
      // expected
    }
  }
}
