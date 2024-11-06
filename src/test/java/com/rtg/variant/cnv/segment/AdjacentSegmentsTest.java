/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.variant.cnv.segment;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class AdjacentSegmentsTest {

  @Test
  public void testConstructor() {
    final Segment segment = new Segment("test", 0, 10, 10.1, 0.4, 1, 1);
    final Segment segment2 = new Segment("test", 10, 40, 10.1, 0.4, 1, 1);
    final AdjacentSegments adjacentSegments = new AdjacentSegments(19.0, segment, segment2);
    assertEquals(segment.getStart(), adjacentSegments.getFirst().getStart());
    assertEquals(segment.getEnd(), adjacentSegments.getFirst().getEnd());

    assertEquals(segment2.getStart(), adjacentSegments.getSecond().getStart());
    assertEquals(segment2.getEnd(), adjacentSegments.getSecond().getEnd());

    assertEquals(19.0, adjacentSegments.getScore(), 1e-8);


  }

}
