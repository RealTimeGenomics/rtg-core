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
    for (long i = 1; i <= mutableGraph.numberPaths(); i++) {
      assertEquals(deleted.contains(i), mutableGraph.pathDeleted(i));
    }
  }

  private MutableGraph baseGraph() {
    final HashMap<String, String> pathAttr = new HashMap<>();
    pathAttr.put(GraphKmerAttribute.READ_COUNT, "count em all");
    MutableGraph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "AAAA", "AAA", "A"}, new long[][]{{1, 2}, {1, 3}, {2, 3}, {3, 4}, {1, 4}}, Collections.<String, String>emptyMap(), pathAttr);
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
    MutableGraph g = simpleCross();
    List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);

    actual = FilterPaths.unambiguousPath(3, 2, g, 10);
    assertEquals(Arrays.asList(3L, 2L, 5L), actual);

    actual = FilterPaths.unambiguousPath(-4, -2, g, 10);
    assertEquals(Arrays.asList(-4L, -2L, -1L), actual);
  }
  public void testUnambiguousImprove() {
    MutableGraph g = simpleCross();
    FilterPaths.improveMultiple(g, 10);
    assertTrue(g.pathDeleted(1));
    assertTrue(g.pathDeleted(2));
    assertTrue(g.pathDeleted(3));
    assertTrue(g.pathDeleted(4));
  }
  public void testUnambiguousExtraPath() {
    MutableGraph g = simpleCross();
    long contigId = g.addContig(new ContigString("AA"));
    g.addPath(new PathArray(1, contigId));
    List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);
  }
  public void testUnambiguousExtraPath2() {
    MutableGraph g = simpleCross();
    long contigId = g.addContig(new ContigString("AA"));
    g.addPath(new PathArray(contigId, 1));
    List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);
  }
  public void testSimplifyExtraPath() {
    MutableGraph g = simpleCross();
    long contigId = g.addContig(new ContigString("AA"));
    long pathId = g.addPath(new PathArray(contigId, 1));
    long longerPath = g.addPath(new PathArray(contigId, 1, 2, 4));
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
    MutableGraph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "AAAA", "AAA", "A", "AA"}, new long[][]{{1, 2}, {3, 2}, {2, 4}, {2, 5}, {1, 2, 4}, {3, 2, 5}}, Collections.<String, String>emptyMap(), pathAttr);
    g.setPathAttribute(5, GraphKmerAttribute.READ_COUNT, "10");
    g.setPathAttribute(6, GraphKmerAttribute.READ_COUNT, "10");
    return g;
  }

  public void testAmbiguousPathSimplest() {
    MutableGraph g = ambiguousCross();
    List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L), actual);

    actual = FilterPaths.unambiguousPath(3, 2, g, 10);
    assertEquals(Arrays.asList(3L, 2L, 5L), actual);
  }
  public void testBelowCoverage() {
    MutableGraph g = ambiguousCross();
    g.setPathAttribute(5, GraphKmerAttribute.READ_COUNT, "9");
    g.setPathAttribute(6, GraphKmerAttribute.READ_COUNT, "9");
    List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(0, actual.size());

  }

  private MutableGraph ambiguousCross() {
    final HashMap<String, String> pathAttr = new HashMap<>();
    pathAttr.put(GraphKmerAttribute.READ_COUNT, "count em all");
    MutableGraph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "AAAA", "AAA", "A", "AA"}, new long[][]{{1, 2}, {3, 2}, {2, 4}, {2, 5}, {1, 2, 4}, {3, 2, 5}, {1, 2, 5}}, Collections.<String, String>emptyMap(), pathAttr);
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
    MutableGraph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "AAAA", "AAA", "A", "AA", "AA"}
        , new long[][]{{1, 2}, {3, 2}, {2, 4}, {4, 5}, {4, 6}, {1, 2, 4}, {3, 2, 4}, {1, 2, 4, 5}, {3, 2, 4, 6}, {3, 2, 4, 5}}, Collections.<String, String>emptyMap(), pathAttr);
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
    MutableGraph g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "AAAA", "AAA", "A", "AA", "AA"}
        , new long[][]{{1, 2}, {3, 2}, {2, 4}, {4, 5}, {4, 6}, {1, 2, 4}, {3, 2, 4}, {1, 2, 4, 5}, {3, 2, 4, 6}, {3, 2, 4, 5}}, Collections.<String, String>emptyMap(), pathAttr);
    g.setPathAttribute(6, GraphKmerAttribute.READ_COUNT, "1");
    g.setPathAttribute(7, GraphKmerAttribute.READ_COUNT, "10");
    g.setPathAttribute(8, GraphKmerAttribute.READ_COUNT, "9");
    g.setPathAttribute(9, GraphKmerAttribute.READ_COUNT, "10");
    List<Long> actual = FilterPaths.unambiguousPath(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);
  }

  public void testBiDirectional() {
    MutableGraph g = simpleCross();
    List<Long> actual = FilterPaths.biDirectional(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);
  }

  public void testBiReverseAmbiguous() {
    MutableGraph g = simpleCross();
    g.addPath(new PathArray(3, 2, 4));
    List<Long> actual = FilterPaths.biDirectional(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L), actual);
  }

  public void testBiDirectionalLonger() {
    MutableGraph g = simpleCross();
    long contigId = g.addContig(new ContigString("GGG"));
    long pathId = g.addPath(new PathArray(contigId, 1, 2, 4));
    g.setPathAttribute(pathId, GraphKmerAttribute.READ_COUNT, "10");
    List<Long> actual = FilterPaths.biDirectional(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L, 4L), actual);
  }
  public void testBiDirectionalAmbiguous() {
    MutableGraph g = ambiguousCross();
    List<Long> actual = FilterPaths.biDirectional(1, 2, g, 10);
    assertEquals(Arrays.asList(1L, 2L), actual);
    actual = FilterPaths.biDirectional(3, 2, g, 10);
    assertEquals(Arrays.asList(3L, 2L), actual);
  }
}
