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
package com.rtg.variant.bayes.multisample.cancer;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
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
