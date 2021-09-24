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

import junit.framework.TestCase;

/**
 */
public class DescriptionCountsTest extends TestCase {

  public void testStuff() {
    final DescriptionCounts dc = new DescriptionCounts(4, 1);
    assertEquals(4, dc.getSize());
    assertEquals(1, dc.getReferenceIndex());

    dc.increment(0, 2);
    dc.increment(1, 1);
    dc.increment(2, 4);
    dc.increment(3, 5);

    assertEquals(12, dc.getTotalCount());

    assertNotNull(dc.getCounts());

    assertEquals(2, dc.getCount(0));
    assertEquals(1, dc.getCount(1));
    assertEquals(4, dc.getCount(2));
    assertEquals(5, dc.getCount(3));
  }
}
