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

package com.rtg.index;

import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Test
 */
public class RepeatFrequencyFilterMethodTest extends TestCase {

  @Override
  protected void setUp() throws Exception {
    super.setUp();
    Diagnostic.setLogStream();
  }

  public void testSimple() {
    final RepeatFrequencyFilterMethod m = new RepeatFrequencyFilterMethod(50, false, 100, 0);
    m.initialize(null);
    assertTrue(m.keepHash(0, 5));
    assertTrue(m.keepHash(0, 25));
    assertTrue(m.keepHash(0, 50));
    assertFalse(m.keepHash(0, 80));
    assertFalse(m.keepHash(0, 1000));
    assertFalse(m.keepHash(0, 10000));
  }

  public void testProportional() {
    final RepeatFrequencyFilterMethod m = new RepeatFrequencyFilterMethod(50, true, 100, 0);
    //discard half the hashes
    final SparseFrequencyHistogram sph = new SparseFrequencyHistogram();
    for (int i = 1; i <= 10; i++) {
      sph.add(i, 10 / i);
    }
    m.internalInitializeProportional(sph, 100);
    assertTrue(m.keepHash(0, 2));
    assertTrue(m.keepHash(0, 3));
    assertTrue(m.keepHash(0, 4));
    assertFalse(m.keepHash(0, 5));
    assertFalse(m.keepHash(0, 6));
    assertFalse(m.keepHash(0, 9));
  }
}