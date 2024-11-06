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

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;

import junit.framework.TestCase;

/**
 */
public class FilterPathsTest extends TestCase {
  public void testFilterSingle() {
    checkThreshold(Arrays.asList(2L), 1);
    checkThreshold(Arrays.asList(2L, 5L), 2);
    checkThreshold(Arrays.asList(2L, 5L, 3L), 3);
  }

  private void checkThreshold(List<Long> deleted, int threshold) {
    final MutableGraph mutableGraph = baseGraph();
    FilterPaths.improveSingle(mutableGraph, threshold);
    for (long i = 1; i <= mutableGraph.numberPaths(); ++i) {
      assertEquals(deleted.contains(i), mutableGraph.pathDeleted(i));
    }
  }

  private MutableGraph baseGraph() {
    final HashMap<String, String> pathAttr = new HashMap<>();
    pathAttr.put(GraphKmerAttribute.READ_COUNT, "count em all");
    final MutableGraph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "AAAA", "AAA", "A"}, new long[][]{{1, 2}, {1, 3}, {2, 3}, {3, 4}, {1, 4}}, Collections.emptyMap(), pathAttr);
    g.setPathAttribute(1, GraphKmerAttribute.READ_COUNT, "10");
    g.setPathAttribute(3, GraphKmerAttribute.READ_COUNT, "2");
    g.setPathAttribute(4, GraphKmerAttribute.READ_COUNT, "5");
    g.setPathAttribute(5, GraphKmerAttribute.READ_COUNT, "1");
    return g;
  }

  // 1   4
  //  \ /
  //   2
  //  / \
  // 3   5
  public void testUnambiguousPathSimplest() {
    final MutableGraph g = simpleCross();
    List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);

    actual = FilterPaths.unambiguousPath(3, 2, g, 10);
    assertEquals(Arrays.asList(3L, 2L, 5L), actual);

    actual = FilterPaths.unambiguousPath(-4, -2, g, 10);
    assertEquals(Arrays.asList(-4L, -2L, -1L), actual);
  }
  public void testUnambiguousImprove() {
    final MutableGraph g = simpleCross();
    FilterPaths.improveMultiple(g, 10);
    assertTrue(g.pathDeleted(1));
    assertTrue(g.pathDeleted(2));
    assertTrue(g.pathDeleted(3));
    assertTrue(g.pathDeleted(4));
  }
  public void testUnambiguousExtraPath() {
    final MutableGraph g = simpleCross();
    final long contigId = g.addContig(new ContigString("AA"));
    g.addPath(new PathArray(1, contigId));
    final List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);
  }
  public void testUnambiguousExtraPath2() {
    final MutableGraph g = simpleCross();
    final long contigId = g.addContig(new ContigString("AA"));
    g.addPath(new PathArray(contigId, 1));
    final List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);
  }
  public void testSimplifyExtraPath() {
    final MutableGraph g = simpleCross();
    final long contigId = g.addContig(new ContigString("AA"));
    final long pathId = g.addPath(new PathArray(contigId, 1));
    final long longerPath = g.addPath(new PathArray(contigId, 1, 2, 4));
    g.setPathAttribute(longerPath, GraphKmerAttribute.READ_COUNT, "5");
    FilterPaths.improveMultiple(g, 5);
    assertTrue(g.pathDeleted(1));
    assertTrue(g.pathDeleted(2));
    assertTrue(g.pathDeleted(3));
    assertTrue(g.pathDeleted(4));
    assertTrue(g.pathDeleted(pathId));
  }

  private MutableGraph simpleCross() {
    final HashMap<String, String> pathAttr = new HashMap<>();
    pathAttr.put(GraphKmerAttribute.READ_COUNT, "count em all");
    final MutableGraph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "AAAA", "AAA", "A", "AA"}, new long[][]{{1, 2}, {3, 2}, {2, 4}, {2, 5}, {1, 2, 4}, {3, 2, 5}}, Collections.emptyMap(), pathAttr);
    g.setPathAttribute(5, GraphKmerAttribute.READ_COUNT, "10");
    g.setPathAttribute(6, GraphKmerAttribute.READ_COUNT, "10");
    return g;
  }

  public void testAmbiguousPathSimplest() {
    final MutableGraph g = ambiguousCross();
    List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L), actual);

    actual = FilterPaths.unambiguousPath(3, 2, g, 10);
    assertEquals(Arrays.asList(3L, 2L, 5L), actual);
  }
  public void testBelowCoverage() {
    final MutableGraph g = ambiguousCross();
    g.setPathAttribute(5, GraphKmerAttribute.READ_COUNT, "9");
    g.setPathAttribute(6, GraphKmerAttribute.READ_COUNT, "9");
    final List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(0, actual.size());

  }

  private MutableGraph ambiguousCross() {
    final HashMap<String, String> pathAttr = new HashMap<>();
    pathAttr.put(GraphKmerAttribute.READ_COUNT, "count em all");
    final MutableGraph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "AAAA", "AAA", "A", "AA"}, new long[][]{{1, 2}, {3, 2}, {2, 4}, {2, 5}, {1, 2, 4}, {3, 2, 5}, {1, 2, 5}}, Collections.emptyMap(), pathAttr);
    g.setPathAttribute(5, GraphKmerAttribute.READ_COUNT, "10");
    g.setPathAttribute(6, GraphKmerAttribute.READ_COUNT, "10");
    return g;
  }

  // 1     5
  //  \   /
  //   2-4
  //  /   \
  // 3     6
  public void testUnambiguousPathLonger() {
    final HashMap<String, String> pathAttr = new HashMap<>();
    pathAttr.put(GraphKmerAttribute.READ_COUNT, "count em all");
    final MutableGraph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "AAAA", "AAA", "A", "AA", "AA"}
        , new long[][]{{1, 2}, {3, 2}, {2, 4}, {4, 5}, {4, 6}, {1, 2, 4}, {3, 2, 4}, {1, 2, 4, 5}, {3, 2, 4, 6}, {3, 2, 4, 5}}, Collections.emptyMap(), pathAttr);
    g.setPathAttribute(6, GraphKmerAttribute.READ_COUNT, "10");
    g.setPathAttribute(7, GraphKmerAttribute.READ_COUNT, "10");
    g.setPathAttribute(8, GraphKmerAttribute.READ_COUNT, "10");
    g.setPathAttribute(9, GraphKmerAttribute.READ_COUNT, "10");
    List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L, 5L), actual);

    actual = FilterPaths.unambiguousPath(3, 2, g, 10);
    assertEquals(Arrays.asList(3L, 2L, 4L), actual);

    actual = FilterPaths.unambiguousPath(-6, -4, g, 10);
    assertEquals(Arrays.asList(-6L, -4L, -2L, -3L), actual);

  }
  public void testUnambiguousPathReadCount() {
    final HashMap<String, String> pathAttr = new HashMap<>();
    pathAttr.put(GraphKmerAttribute.READ_COUNT, "count em all");
    final MutableGraph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "AAAA", "AAA", "A", "AA", "AA"}
        , new long[][]{{1, 2}, {3, 2}, {2, 4}, {4, 5}, {4, 6}, {1, 2, 4}, {3, 2, 4}, {1, 2, 4, 5}, {3, 2, 4, 6}, {3, 2, 4, 5}}, Collections.emptyMap(), pathAttr);
    g.setPathAttribute(6, GraphKmerAttribute.READ_COUNT, "1");
    g.setPathAttribute(7, GraphKmerAttribute.READ_COUNT, "10");
    g.setPathAttribute(8, GraphKmerAttribute.READ_COUNT, "9");
    g.setPathAttribute(9, GraphKmerAttribute.READ_COUNT, "10");
    final List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);
  }

  public void testBiDirectional() {
    final MutableGraph g = simpleCross();
    final List<Long> actual = FilterPaths.biDirectional(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);
  }

  public void testBiReverseAmbiguous() {
    final MutableGraph g = simpleCross();
    g.addPath(new PathArray(3, 2, 4));
    final List<Long> actual = FilterPaths.biDirectional(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L), actual);
  }

  public void testBiDirectionalLonger() {
    final MutableGraph g = simpleCross();
    final long contigId = g.addContig(new ContigString("GGG"));
    final long pathId = g.addPath(new PathArray(contigId, 1, 2, 4));
    g.setPathAttribute(pathId, GraphKmerAttribute.READ_COUNT, "10");
    final List<Long> actual = FilterPaths.biDirectional(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);
  }
  public void testBiDirectionalAmbiguous() {
    final MutableGraph g = ambiguousCross();
    List<Long> actual = FilterPaths.biDirectional(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L), actual);
    actual = FilterPaths.biDirectional(3, 2, g, 10);
    assertEquals(Arrays.asList(3L, 2L), actual);
  }
}
