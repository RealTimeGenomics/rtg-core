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

package com.rtg.variant.bayes.multisample.population;

import java.util.HashMap;
import java.util.Set;

import junit.framework.TestCase;

/**
 */
public class AlleleCountsTest extends TestCase {

  public void testEmpty() {
    final HashMap<String, Integer> countsMap = new HashMap<>();
    final AlleleCounts ac = new AlleleCounts(0, countsMap, "A");
    assertEquals(0, ac.position());
    assertEquals(1, ac.refLength());
    assertFalse(ac.isComplex());
    assertEquals("A", ac.getReferenceAllele());
    assertEquals(-1, ac.count("AC"));
    assertEquals(1, ac.maxLength());
    assertEquals("(0 {})", ac.toString());
  }

  public void test() {
    final HashMap<String, Integer> countsMap = new HashMap<>();
    countsMap.put("AC", 123);
    countsMap.put("T", 22);
    final AlleleCounts ac = new AlleleCounts(456, countsMap, "AC");
    assertEquals(456, ac.position());
    assertEquals(2, ac.refLength());
    assertTrue(ac.isComplex());
    assertEquals("AC", ac.getReferenceAllele());
    final Set<String> seen = ac.allelesSeen();
    assertNotNull(seen);
    assertEquals(2, seen.size());
    assertEquals(123, ac.count("AC"));
    assertEquals(22, ac.count("T"));
    assertEquals(2, ac.maxLength());
    assertEquals("(456 {AC=123, T=22})", ac.toString());
  }

  public void testNotComplex() {
    final HashMap<String, Integer> countsMap = new HashMap<>();
    countsMap.put("A", 1);
    countsMap.put("T", 1);
    final AlleleCounts ac = new AlleleCounts(456, countsMap, "A");
    assertEquals(1, ac.refLength());
    assertFalse(ac.isComplex());
    assertEquals(1, ac.maxLength());
  }

  public void testRefOnlyComplex() {
    final HashMap<String, Integer> countsMap = new HashMap<>();
    countsMap.put("A", 1);
    countsMap.put("T", 1);
    final AlleleCounts ac = new AlleleCounts(456, countsMap, "AC");
    assertEquals(2, ac.refLength());
    assertTrue(ac.isComplex());
    assertEquals(2, ac.maxLength());
  }

  public void testNonRefOnlyComplex() {
    final HashMap<String, Integer> countsMap = new HashMap<>();
    countsMap.put("AC", 1);
    countsMap.put("T", 1);
    final AlleleCounts ac = new AlleleCounts(456, countsMap, "A");
    assertEquals(1, ac.refLength());
    assertTrue(ac.isComplex());
    assertEquals(2, ac.maxLength());
  }

  public void testZeroLength() {
    final HashMap<String, Integer> countsMap = new HashMap<>();
    countsMap.put("A", 1);
    countsMap.put("", 1);
    final AlleleCounts ac = new AlleleCounts(456, countsMap, "A");
    assertEquals(1, ac.refLength());
    assertTrue(ac.isComplex());
    assertEquals(1, ac.maxLength());
  }
}
