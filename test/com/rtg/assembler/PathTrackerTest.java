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
    List<Long> result = new ArrayList<>();
    for (long l : longs) {
      result.add(l);
    }
    return result;
  }
  public void testTracker() {
    Graph graph = GraphMapCliTest.makeGraph(0, new String[] {"ACCC", "TTGG", "CGGGG", "ACGT", "AACT", "ATTGTTAACAAT", "ACGCGA", "ACCGGT"}, new long[][] {});
    PathTracker tracker = new PathTracker(new PalindromeTracker(graph));
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
    MutableGraph graph = GraphMapCliTest.makeGraph(0, new String[] {"ACCC", "TTGG", "CGGGG"}, new long[][] {{1, 2, 3}});
    graph.addPathAttribute(READ_COUNT, GraphKmerAttribute.READ_COUNT_DESCRIPTION);
    graph.setPathAttribute(1, READ_COUNT, "100");
    PathTracker tracker1 = new PathTracker(new PalindromeTracker(graph));
    PathTracker tracker2 = new PathTracker(new PalindromeTracker(graph));

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
