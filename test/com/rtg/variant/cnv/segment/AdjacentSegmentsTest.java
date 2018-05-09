/*
 * Copyright (c) 2016. Real Time Genomics Limited.
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

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class AdjacentSegmentsTest {

  @Test
  public void testConstructor() {
    final Segment segment = new Segment("test", 0, 10, 10.1, 0.4);
    final Segment segment2 = new Segment("test", 10, 40, 10.1, 0.4);
    final AdjacentSegments adjacentSegments = new AdjacentSegments(19.0, segment, segment2);
    assertEquals(segment.getStart(), adjacentSegments.getFirst().getStart());
    assertEquals(segment.getEnd(), adjacentSegments.getFirst().getEnd());

    assertEquals(segment2.getStart(), adjacentSegments.getSecond().getStart());
    assertEquals(segment2.getEnd(), adjacentSegments.getSecond().getEnd());

    assertEquals(19.0, adjacentSegments.getScore(), 1e-8);


  }

}
