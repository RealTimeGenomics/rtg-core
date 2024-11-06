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

package com.rtg.assembler;

import static com.rtg.assembler.graph.implementation.GraphKmerAttribute.READ_COUNT;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;

import junit.framework.TestCase;

/**
 */
public class PathTrackerTest extends TestCase {
  List<Long> longs(long ... longs) {
    final List<Long> result = new ArrayList<>();
    for (long l : longs) {
      result.add(l);
    }
    return result;
  }
  public void testTracker() {
    final Graph graph = GraphMapCliTest.makeGraph(0, new String[] {"ACCC", "TTGG", "CGGGG", "ACGT", "AACT", "ATTGTTAACAAT", "ACGCGA", "ACCGGT"}, new long[][] {});
    final PathTracker tracker = new PathTracker(new PalindromeTracker(graph));
    tracker.increment(longs(1L, 2L, 3L));
    tracker.increment(longs(-3L, -2L, -1L));
    tracker.increment(longs(1L, 2L, 3L));

    assertEquals(1, tracker.mPathCounts.entrySet().size());
    for (Map.Entry<List<Long>, Integer> entry : tracker.mPathCounts.entrySet()) {
      assertEquals(Arrays.asList(-3L, -2L, -1L), entry.getKey());
      assertEquals(Integer.valueOf(3), entry.getValue());
    }

    tracker.increment(longs(1, 4, 3));
    tracker.increment(longs(1, -4, 3));
    tracker.increment(longs(-3, 4, -1));
    tracker.increment(longs(-3, -4, -1));
    assertEquals(Integer.valueOf(4), tracker.mPathCounts.get(longs(-3, 4, -1)));

    tracker.increment(longs(-3, -4, 1));
    assertEquals(Integer.valueOf(4), tracker.mPathCounts.get(longs(-3, 4, -1)));
    assertEquals(Integer.valueOf(1), tracker.mPathCounts.get(longs(-3, 4, 1)));

    tracker.increment(longs(4, 6, 8));
    tracker.increment(longs(-4, -6, -8));
    tracker.increment(longs(-8, -6, -4));
    tracker.increment(longs(8, 6, 4));
    tracker.increment(longs(8, -6, 4));
    tracker.increment(longs(-8, -6, 4));

    assertEquals(Integer.valueOf(6), tracker.mPathCounts.get(longs(4, 6, 8)));

    tracker.increment(longs(4, 4, 4));
    tracker.increment(longs(4, 4, -4));
    assertEquals(Integer.valueOf(2), tracker.mPathCounts.get(longs(4, 4, 4)));
  }

  public void testMerge() {
    final MutableGraph graph = GraphMapCliTest.makeGraph(0, new String[] {"ACCC", "TTGG", "CGGGG"}, new long[][] {{1, 2, 3}});
    graph.addPathAttribute(READ_COUNT, GraphKmerAttribute.READ_COUNT_DESCRIPTION);
    graph.setPathAttribute(1, READ_COUNT, "100");
    final PathTracker tracker1 = new PathTracker(new PalindromeTracker(graph));
    final PathTracker tracker2 = new PathTracker(new PalindromeTracker(graph));

    tracker1.increment(longs(1L, 2L, 3L));
    tracker1.increment(longs(1L, 3L, 3L));

    tracker2.increment(longs(1L, 2L, 3L));
    tracker2.increment(longs(1L, 2L, 2L));

    final SortedMap<List<Long>, Integer> merged = PathTracker.merge(Arrays.asList(tracker1, tracker2));

    assertEquals(Integer.valueOf(2), merged.get(longs(-3, -2, -1)));

    assertEquals(Integer.valueOf(1), merged.get(longs(-2, -2, -1)));
    assertEquals(Integer.valueOf(1), merged.get(longs(-3, -3, -1)));

    PathTracker.apply(merged, graph);
    assertEquals("102", graph.pathAttribute(1, READ_COUNT));

    assertEquals(-3, graph.pathContig(2, 0));
    assertEquals(-3, graph.pathContig(2, 1));
    assertEquals(-1, graph.pathContig(2, 2));

    assertEquals(-2, graph.pathContig(3, 0));
    assertEquals(-2, graph.pathContig(3, 1));
    assertEquals(-1, graph.pathContig(3, 2));

    assertEquals("1", graph.pathAttribute(2, READ_COUNT));
    assertEquals("1", graph.pathAttribute(3, READ_COUNT));
  }


}
