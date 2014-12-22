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
