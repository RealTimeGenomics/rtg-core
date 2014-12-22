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

import java.util.HashSet;
import java.util.Set;

import com.rtg.assembler.graph.Graph;

import junit.framework.TestCase;

/**
 */
public class GraphTraversionsTest extends TestCase {
  Set<Long> longs(long... list) {
    Set<Long> result = new HashSet<>();
    for (long l :list) {
      result.add(l);
    }
    return result;
  }
  public void test() {
    Graph g = GraphMapCliTest.makeGraph(0, new String[]{"ACGT", "AACC", "ACAA", "GGAT", "ACGGT"}, new long[][]{{1, 2, 3}, {1, 2}, {-4, -3}, {3, 5}});
    final GraphTraversions traversions = new GraphTraversions(g);
    assertEquals(longs(4, 5), traversions.get(3L).next());
    assertEquals(longs(2), traversions.get(3L).previous());
    assertEquals(longs(-4, -5), traversions.get(-3L).previous());
  }

  public void testPalindrome() {
    Graph g = GraphMapCliTest.makeGraph(0, new String[]{"ACGT", "AACC", "ACAA", "GGCC", "ACGGT"}, new long[][]{{-1, 2}, {1, 2}, {3, 4}, {3, -4}, {4, 5}, {-4, 5}});
    final GraphTraversions traversions = new GraphTraversions(g);
    assertEquals(longs(1, -1), traversions.get(2L).previous());
    assertEquals(longs(1, -1), traversions.get(-2L).next());

    assertEquals(longs(4, -4), traversions.get(3L).next());
    assertEquals(longs(4, -4), traversions.get(5L).previous());
    assertEquals(longs(-3, 5), traversions.get(4L).next());
    assertEquals(longs(3, -5), traversions.get(4L).previous());
  }

}
