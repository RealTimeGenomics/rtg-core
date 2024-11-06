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

import java.util.HashSet;
import java.util.Set;

import com.rtg.assembler.graph.Graph;

import junit.framework.TestCase;

/**
 */
public class GraphTraversionsTest extends TestCase {
  Set<Long> longs(long... list) {
    final Set<Long> result = new HashSet<>();
    for (long l :list) {
      result.add(l);
    }
    return result;
  }
  public void test() {
    final Graph g = GraphMapCliTest.makeGraph(0, new String[]{"ACGT", "AACC", "ACAA", "GGAT", "ACGGT"}, new long[][]{{1, 2, 3}, {1, 2}, {-4, -3}, {3, 5}});
    final GraphTraversions traversions = new GraphTraversions(g);
    assertEquals(longs(4, 5), traversions.get(3L).next());
    assertEquals(longs(2), traversions.get(3L).previous());
    assertEquals(longs(-4, -5), traversions.get(-3L).previous());
  }

  public void testPalindrome() {
    final Graph g = GraphMapCliTest.makeGraph(0, new String[]{"ACGT", "AACC", "ACAA", "GGCC", "ACGGT"}, new long[][]{{-1, 2}, {1, 2}, {3, 4}, {3, -4}, {4, 5}, {-4, 5}});
    final GraphTraversions traversions = new GraphTraversions(g);
    assertEquals(longs(1, -1), traversions.get(2L).previous());
    assertEquals(longs(1, -1), traversions.get(-2L).next());

    assertEquals(longs(4, -4), traversions.get(3L).next());
    assertEquals(longs(4, -4), traversions.get(5L).previous());
    assertEquals(longs(-3, 5), traversions.get(4L).next());
    assertEquals(longs(3, -5), traversions.get(4L).previous());
  }

}
