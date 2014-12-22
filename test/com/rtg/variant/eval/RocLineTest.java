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

package com.rtg.variant.eval;

import junit.framework.TestCase;

/**
 */
public class RocLineTest extends TestCase {


  void check(RocLine roc1, RocLine roc2) {
    assertTrue(roc1.compareTo(roc2) < 0);
    assertTrue(roc2.compareTo(roc1) > 0);
    assertTrue(roc1.compareTo(roc1) == 0);
    assertTrue(roc2.compareTo(roc2) == 0);
  }
  public void testPosterior() {
    final RocLine roc1 = new RocLine("1", 2, 3.5, 4, false);
    final RocLine roc2 = new RocLine("1", 2, 3, 4, false);
    check(roc1, roc2);
  }
  public void testSequence() {
    final RocLine roc1 = new RocLine("1", 2, 3, 4, false);
    final RocLine roc2 = new RocLine("2", 2, 3, 4, false);
    check(roc1, roc2);
  }
  public void testPosition() {
    final RocLine roc1 = new RocLine("1", 2, 3, 4, false);
    final RocLine roc2 = new RocLine("1", 3, 3, 4, false);
    check(roc1, roc2);
  }

  public void testCorrect() {
    final RocLine roc1 = new RocLine("1", 2, 3, 4, false);
    final RocLine roc2 = new RocLine("1", 2, 3, 4, true);
    check(roc1, roc2);
  }

  public void testFieldOrder() {
    RocLine roc1 = new RocLine("3", 6, 2, 4, true);
    final RocLine roc2 = new RocLine("2", 2, 6, 4, false);
    check(roc2, roc1);
    roc1 = new RocLine("3", 6, 7, 4, true);
    check(roc1, roc2);
    roc1 = new RocLine("1", 6, 6, 4, false);
    check(roc1, roc2);
    roc1 = new RocLine("1", 1, 6, 4, false);
    check(roc1, roc2);
  }
  public void testEquals() {
    final RocLine roc1 = new RocLine("3", 6, 2, 4, true);
    final RocLine roc2 = new RocLine("2", 2, 6, 4, false);
    assertFalse(roc1.equals(null));
    assertFalse(roc1.equals(roc2));
    assertTrue(roc1.equals(roc1));
  }

}
