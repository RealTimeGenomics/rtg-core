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
package com.rtg.variant;

import junit.framework.TestCase;

/**
 *         Date: 27/10/11
 *         Time: 4:39 PM
 */
public class StaticThresholdTest extends TestCase {
  public void testThreshold() {
        final StaticThreshold threshold = new StaticThreshold(60);
        assertEquals(60, threshold.thresholdSingle("chr1"));
        assertEquals(60, threshold.thresholdTotal("chr2"));
        assertEquals("60:60", threshold.toString());

  }
}
