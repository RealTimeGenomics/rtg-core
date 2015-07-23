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
package com.rtg.simulation;

import junit.framework.TestCase;

/**
 */
public class ReadSimEvalStatisticsTest extends TestCase {


  public void test() {
    final ReadSimEvalStatistics s = new ReadSimEvalStatistics(20);
    assertEquals(s.length(), 20);
    s.found(0);
    assertTrue(s.isFound(0));
    s.mated(0);

    assertTrue(s.isMated(0));

    s.better(1);
    assertTrue(s.isBetter(1));
    assertFalse(s.isFound(1));
    //assertFalse(s.isBetter(1));
    assertFalse(s.isMultiple(1));
    assertFalse(s.isUnmapped(1));
    assertFalse(s.isUnmated(1));

    s.multiple(1);
    assertTrue(s.isMultiple(1));
    s.unmapped(2);
    assertTrue(s.isUnmapped(2));
    s.unmated(3);
    assertTrue(s.isUnmated(3));

    assertFalse(s.isMapped(4));
    s.mapped(4);
    assertTrue(s.isMapped(4));

  }

}
