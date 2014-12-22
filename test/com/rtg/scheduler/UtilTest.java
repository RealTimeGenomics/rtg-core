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

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import com.rtg.util.integrity.Exam.ExamException;

import junit.framework.TestCase;

/**
 */
public class UtilTest extends TestCase {

  public void testCheckOrder() {
    final Set<Integer> set = new HashSet<>();
    set.add(4);
    set.add(20);
    set.add(null);
    checkCheckOrder(3, set, +1);
    checkCheckOrder(23, set, -1);
  }

  public void testCheckOrderEmpty() {
    final Set<Integer> set = new HashSet<>();
    Util.checkOrder(3, set, +1);
    Util.checkOrder(3, set, -1);
  }

  <X> void checkCheckOrder(final Comparable<X> x0, final Collection<X> set, final int exp) {
    assertTrue(Util.checkOrder(x0, set, exp));
    try {
      Util.checkOrder(x0, set, -exp);
      fail();
    } catch (final ExamException e) {
      // expected
    }
  }
  public void testNonNullSize() {
    final Set<Integer> set = new HashSet<>();
    set.add(4);
    set.add(20);
    set.add(null);
    assertEquals(2, Util.nonNullSize(set));
    assertEquals(0, Util.nonNullSize(new HashSet<Integer>()));
  }

}
