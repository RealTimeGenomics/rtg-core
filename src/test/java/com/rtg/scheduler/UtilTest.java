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
