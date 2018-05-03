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

package com.rtg.variant.cnv.segment;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class EnergySegmentScorerTest extends TestCase {

  public void test() {
    final Segment s = new Segment(0, 1, 100, 1);
    final Segment t = new Segment(new Segment(new Segment(1, 2, 50, 1), new Segment(2, 3, 50, 10)), new Segment(3, 4, 50, 20));
    assertEquals(1875, new EnergySegmentScorer(0, 0).score(s, t), 1e-4);
    assertEquals(1877.772588, new EnergySegmentScorer(1, 0).score(s, t), 1e-4);
    assertEquals(1878.46574, new EnergySegmentScorer(1, 1).score(s, t), 1e-4);
  }
}
