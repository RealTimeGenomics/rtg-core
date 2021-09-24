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

import java.util.Arrays;

import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Test
 */
public class BlacklistFilterMethodTest extends TestCase {

  @Override
  protected void setUp() throws Exception {
    super.setUp();
    Diagnostic.setLogStream();
  }

  public void test() {
    final BlacklistFilterMethod m = new BlacklistFilterMethod(Arrays.asList(0b1111L, 0b1010L, 0b1100L), 4, 1);
    m.initialize(null);
    assertFalse(m.keepHash(0b1111L, 10000));
    assertTrue(m.keepHash(0b1110L, 10000));
    assertFalse(m.keepHash(0b1100L, 10000));
    assertTrue(m.keepHash(0b1101L, 10000));
    assertFalse(m.keepHash(0b1010L, 10000));
    assertTrue(m.keepHash(0b1011L, 10000));
  }
}