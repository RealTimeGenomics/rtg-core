package com.rtg.variant.bayes.multisample.cancer;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 *
 * @author Sean A. Irvine
 */
public class TotalScorerTest extends TestCase {

  public void test() {
    final Scorer ts = new TotalScorer();
    ts.add(null, 42, 42);
    assertEquals(0, ts.size());
    ts.add(2.0, 2, 3);
    assertEquals(1, ts.size());
    ts.add(1.0, 1, 2);
    ts.add(4.0, 4, 5);
    assertEquals(3, ts.size());
    ts.add(null, 42, 45);
    ts.add(1.0, 56, 57);
    assertEquals(4, ts.size());
    assertEquals(63, ts.getTotalRefCount());
    assertEquals(67, ts.getTotalAltCount());
  }

}
