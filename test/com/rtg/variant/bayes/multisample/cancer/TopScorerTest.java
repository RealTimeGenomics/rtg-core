package com.rtg.variant.bayes.multisample.cancer;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 *
 * @author Sean A. Irvine
 */
public class TopScorerTest extends TestCase {

  public void test() {
    final TopScorer ts = new TopScorer(3);
    ts.add(null, 42, 42);
    assertEquals(0, ts.size());
    ts.add(2.0, 2, 3);
    assertEquals(1, ts.size());
    ts.add(1.0, 1, 2);
    ts.add(4.0, 4, 5);
    assertEquals(3, ts.size());
    ts.add(null, 42, 45);
    ts.add(1.0, 56, 57);
    assertEquals(3, ts.size());
    assertEquals(4, ts.getRefCount(0));
    assertEquals(5, ts.getAltCount(0));
    assertEquals(2, ts.getRefCount(1));
    assertEquals(3, ts.getAltCount(1));
    assertEquals(1, ts.getRefCount(2));
    assertEquals(2, ts.getAltCount(2));
  }

}
