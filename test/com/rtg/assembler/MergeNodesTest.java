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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.mode.DnaUtils;
import com.rtg.util.store.StoreDirString;

import junit.framework.TestCase;

/**
 */
public class MergeNodesTest extends TestCase {
  List<Long> longs(long ... longs) {
    final ArrayList<Long> result = new ArrayList<>();
    for (long l : longs) {
      result.add(l);
    }
    return result;
  }
  public void testNewPath() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(2, new String[]{"A", "A", "A", "A", "A", "A", "A", "A"}, new long[][]{});
    assertEquals(longs(1, -77, 4), MergeNodes.newPath(-77, longs(2, 3), longs(1, 2, 3, 4), graph));

    assertEquals(longs(1, -77), MergeNodes.newPath(-77, longs(2, 3, 4),          longs(1, 2, 3, 4), graph));
    assertEquals(longs(1, -77), MergeNodes.newPath(-77, longs(2, 3, 4, 5),       longs(1, 2, 3, 4), graph));
    assertEquals(longs(-77, 4), MergeNodes.newPath(-77, longs(1, 2, 3),          longs(1, 2, 3, 4), graph));
    assertEquals(longs(-77, 4), MergeNodes.newPath(-77, longs(6, 1, 2, 3),       longs(1, 2, 3, 4), graph));
    assertEquals(longs(-77),    MergeNodes.newPath(-77, longs(1, 2, 3, 4),       longs(1, 2, 3, 4), graph));
    assertEquals(longs(-77),    MergeNodes.newPath(-77, longs(6, 1, 2, 3, 4, 5), longs(1, 2, 3, 4), graph));
    assertEquals(longs(6, -77, -77, 4), MergeNodes.newPath(-77, longs(1, 2, 3),  longs(6, 1, 2, 3, 1, 2, 3, 4), graph));
    assertEquals(longs(-77, 3), MergeNodes.newPath(-77, longs(1, 1, 2), longs(1, 2, 3), graph));

    assertEquals(longs(9, 2, 3), MergeNodes.newPath(-77, longs(1, 2, 3), longs(9, 2, 3), graph));
    assertEquals(longs(1, 2, 10), MergeNodes.newPath(-77, longs(1, 2, 3), longs(1, 2, 10), graph));
    assertEquals(longs(9, 2, 10), MergeNodes.newPath(-77, longs(1, 2, 3), longs(9, 2, 10), graph));

    assertEquals(longs(1, 2, 2, 3), MergeNodes.newPath(-77, longs(1, 2, 3), longs(1, 2, 2, 3), graph));
    assertEquals(longs(2, 4), MergeNodes.newPath(-77, longs(1, 2, 3), longs(2, 4), graph));
    assertEquals(longs(4, 2), MergeNodes.newPath(-77, longs(1, 2, 3), longs(4, 2), graph));
  }


  public void testLongestOverlapStart() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(2, new String[]{"A", "A", "A", "A", "A", "A", "A", "A"}, new long[][]{});
    assertEquals(0, MergeNodes.longestStartOverlap(longs(2, 3), longs(1, 2, 3, 4), graph));
    assertEquals(3, MergeNodes.longestStartOverlap(longs(1, 2, 3), longs(1, 2, 3, 4), graph));
    assertEquals(2, MergeNodes.longestStartOverlap(longs(1, 2, 3), longs(2, 3, 4), graph));
    assertEquals(1, MergeNodes.longestStartOverlap(longs(1, 2, 3), longs(3, 4), graph));
    assertEquals(0, MergeNodes.longestStartOverlap(longs(1, 2, 3), longs(4), graph));
    assertEquals(0, MergeNodes.longestStartOverlap(longs(1, 2, 3), longs(1, 2, 4), graph));
    assertEquals(0, MergeNodes.longestStartOverlap(longs(1, 2, 3), longs(2, 2, 3), graph));
    assertEquals(2, MergeNodes.longestStartOverlap(longs(1, 1, 2), longs(1, 2, 3, 4), graph));
    assertEquals(0, MergeNodes.longestStartOverlap(longs(1, 2, 3), longs(6, 1, 2, 3, 4), graph));

    assertEquals(2, MergeNodes.longestStartOverlap(longs(1, 2, 3, 4), longs(2, 3), graph));
  }

  public void testInternalMatch() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(2, new String[]{"A", "A", "A", "A", "A", "A", "A", "A"}, new long[][]{});
    assertTrue(MergeNodes.internalMatch(1, longs(1, 2, 3), longs(6, 1, 2, 3, 4), graph));
    assertFalse(MergeNodes.internalMatch(0, longs(1, 2, 3), longs(6, 1, 2, 3, 4), graph));
    assertFalse(MergeNodes.internalMatch(2, longs(1, 2, 3), longs(6, 1, 2, 3, 4), graph));

    assertTrue(MergeNodes.internalMatch(2, longs(1, 2, 3), longs(6, 7, 1, 2), graph));

    assertFalse(MergeNodes.internalMatch(2, longs(1, 2, 3), longs(6, 7, 1, 2, 4, 5), graph));

  }
  public void testMergeNode() {
    final long[][] paths = {{1, 2, 3, 4, 5}, {2, 3, 4}, {3, 4, 5}, {6, 3, 7}};
    final String[] contigs = {"AAAACC", "CCGGT", "GTTTAT", "ATCCTG", "TGACCAC", "ATAGT", "ATCACAC"};
    final MutableGraph graph = GraphMapCliTest.makeGraph(2, contigs
        , paths);
    final long newId = contigs.length + 1;
    assertEquals(newId, MergeNodes.mergeNodes(graph, 2, longs(2, 3, 4)));
    assertEquals("CCGGTTTATCCTG", ContigString.contigSequenceString(graph.contig(newId)));
    assertTrue(graph.pathDeleted(1));
    assertTrue(graph.pathDeleted(2));
    assertTrue(graph.pathDeleted(3));
    assertFalse(graph.pathDeleted(4));
    final Set<List<Long>> newPaths = new HashSet<>();
    for (int i = paths.length + 1; i <= graph.numberPaths(); ++i) {
      newPaths.add(MergeNodes.getPath(graph, i));
    }
    final Set<List<Long>> expected = new HashSet<>();
    expected.add(longs(1, newId, 5));
    expected.add(longs(newId, 5));
    assertEquals(expected, newPaths);
  }

  void addPath(List<Long> contigs, int readCount, MutableGraph graph) {
    final long id = graph.addPath(new PathArray(contigs));
    graph.setPathAttribute(id, GraphKmerAttribute.READ_COUNT, "" + readCount);
  }

  public void testSimplify() {
    final long[][] paths = {{1, 2}, {2, 3}, {4, 2}, {2, 5}};
    final String[] contigs = {"AAAACC", "CCGGT", "GTTTAT", "ATCC", "GTCCAC"};
    final Map <String, String> attr = new HashMap<>();
    attr.put(GraphKmerAttribute.READ_COUNT, "foo");
    final MutableGraph graph = GraphMapCliTest.makeGraph(2, contigs
        , paths, attr, attr);
    addPath(longs(1, 2, 3), 2, graph);
    addPath(longs(4, 2, 5), 2, graph);
    final MergeNodes merge = new MergeNodes(graph, 2, 2);
    merge.simplifyGraph();
    final Graph compact = graph.compact();

    final String[] expected = {"AAAACCGGTTTAT", "ATCCGGTCCAC"};
    finalContigs(expected, compact);
    assertEquals(0, compact.numberPaths());
  }

  public void testUnMergeableExtraPath() {
    final long[][] paths = {{1, 2}, {2, 3}, {4, 2}, {2, 5}};
    final String[] contigs = {"AAAACC", "CCGGT", "GTTTAT", "ATCC", "GTCCAC"};
    final Map <String, String> attr = new HashMap<>();
    attr.put(GraphKmerAttribute.READ_COUNT, "foo");
    final MutableGraph graph = GraphMapCliTest.makeGraph(2, contigs
        , paths, attr, attr);
    addPath(longs(1, 2, 3), 2, graph);
    addPath(longs(4, 2, 5), 2, graph);
    addPath(longs(4, 2, 3), 2, graph);
    final MergeNodes merge = new MergeNodes(graph, 2, 2);
    merge.simplifyGraph();
    contigDeleted(graph, longs());
    final String[] expected = {};
    newContigs(contigs.length, expected, graph);
    pathDeleted(graph, longs());
  }
  public void testUnMergeableExtraNodes() {
    final long[][] paths = {{1, 2}, {2, 3}, {4, 2}, {2, 5}, {1, 6}, {7, 3}};
    final String[] contigs = {"AAAACC", "CCGGT", "GTTTAT", "ATCC", "GTCCAC", "CCGTG", "AGTGT"};
    final Map <String, String> attr = new HashMap<>();
    attr.put(GraphKmerAttribute.READ_COUNT, "foo");
    final MutableGraph graph = GraphMapCliTest.makeGraph(2, contigs
        , paths, attr, attr);
    addPath(longs(1, 2, 3), 2, graph);
    addPath(longs(4, 2, 5), 2, graph);
    final MergeNodes merge = new MergeNodes(graph, 2, 2);
    merge.simplifyGraph();
    contigDeleted(graph, longs(4, 5));
    final String[] expected = {"ATCCGGTCCAC"};
    newContigs(contigs.length, expected, graph);
    pathDeleted(graph, longs(3, 4, 8));
  }
  public void testUnMergeableExtraNodeLeft() {
    final long[][] paths = {{1, 2}, {2, 3}, {4, 2}, {2, 5}, {1, 6}};
    final String[] contigs = {"AAAACC", "CCGGT", "GTTTAT", "ATCC", "GTCCAC", "CCGTG"};
    final Map <String, String> attr = new HashMap<>();
    attr.put(GraphKmerAttribute.READ_COUNT, "foo");
    final MutableGraph graph = GraphMapCliTest.makeGraph(2, contigs
        , paths, attr, attr);
    addPath(longs(1, 2, 3), 2, graph);
    addPath(longs(4, 2, 5), 2, graph);
    final MergeNodes merge = new MergeNodes(graph, 2, 2);
    merge.simplifyGraph();
    contigDeleted(graph, longs(4, 5));
    final String[] expected = {"ATCCGGTCCAC"};
    newContigs(contigs.length, expected, graph);
    pathDeleted(graph, longs(3, 4, 7));
  }
  public void testUnMergeableExtraNodeRight() {
    final long[][] paths = {{1, 2}, {2, 3}, {4, 2}, {2, 5}, {6, 3}};
    final String[] contigs = {"AAAACC", "CCGGT", "GTTTAT", "ATCC", "GTCCAC", "CCGTG"};
    final Map <String, String> attr = new HashMap<>();
    attr.put(GraphKmerAttribute.READ_COUNT, "foo");
    final MutableGraph graph = GraphMapCliTest.makeGraph(2, contigs
        , paths, attr, attr);
    addPath(longs(1, 2, 3), 2, graph);
    addPath(longs(4, 2, 5), 2, graph);
    final MergeNodes merge = new MergeNodes(graph, 2, 2);
    merge.simplifyGraph();
    contigDeleted(graph, longs(4, 5));
    final String[] expected = {"ATCCGGTCCAC"};
    newContigs(contigs.length, expected, graph);
    pathDeleted(graph, longs(3, 4, 7));
  }

  public void testMergeableExtraNodeRight() {
    final long[][] paths = {{1, 2}, {2, 3}, {4, 2}, {2, 5}, {3, 6}};
    final String[] contigs = {"AAAACC", "CCGGT", "GTTTAT", "ATCC", "GTCCAC", "ATGTG"};
    final Map <String, String> attr = new HashMap<>();
    attr.put(GraphKmerAttribute.READ_COUNT, "foo");
    attr.put(Consensus.COMBINED, "foo");
    attr.put("monkey", "foo");
    final MutableGraph graph = GraphMapCliTest.makeGraph(2, contigs
        , paths, attr, attr);
    graph.setPathAttribute(5, GraphKmerAttribute.READ_COUNT, "4");
    graph.setContigAttribute(1, GraphKmerAttribute.K_MER_FREQ, "2");
    graph.setContigAttribute(2, GraphKmerAttribute.K_MER_FREQ, "4");
    graph.setContigAttribute(1, GraphKmerAttribute.READ_COUNT, "3");
    graph.setContigAttribute(3, GraphKmerAttribute.READ_COUNT, "7");
    graph.setContigAttribute(2, "monkey", "hello");
    addPath(longs(1, 2, 3), 1, graph);
    addPath(longs(1, 2, 3, 6), 1, graph);
    addPath(longs(4, 2, 5), 2, graph);
    final MergeNodes merge = new MergeNodes(graph, 2, 2);
    merge.simplifyGraph();
    final Graph compact = graph.compact();
    final String[] expected = {"AAAACCGGTTTATGTG", "ATCCGGTCCAC"};
    finalContigs(expected, graph);

    assertEquals(0, compact.numberPaths());
    assertEquals("10", compact.contigAttribute(2, GraphKmerAttribute.READ_COUNT));
    assertEquals("-6/-7:(1/2/3)", compact.contigAttribute(2, Consensus.COMBINED));
    assertEquals("hello", compact.contigAttribute(2, "monkey"));

    assertNull(compact.contigAttribute(1, GraphKmerAttribute.READ_COUNT));
    assertEquals("4/8:(2/5)", compact.contigAttribute(1, Consensus.COMBINED));
    assertEquals("hello", compact.contigAttribute(1, "monkey"));
  }
  void finalContigs(String[] contigs, Graph g) {
    final List<String> actual = new ArrayList<>();
    for (long i = 1; i <= g.numberContigs(); ++i) {
      if (!g.contigDeleted(i)) {
        actual.add(canonicalContig(ContigString.contigSequenceString(g.contig(i))));
      }
    }
    assertEquals("Expected <" + contigs.length + "> contigs but there was <" + actual.size() + ">", contigs.length, actual.size());
    final List<String> expected = new ArrayList<>();
    for (String s : contigs) {
      expected.add(canonicalContig(s));
    }
    Collections.sort(expected);
    Collections.sort(actual);
    assertEquals(expected, actual);
  }

  String canonicalContig(String s) {
    final String reverse = DnaUtils.reverseComplement(s);
    if (s.compareTo(reverse) <= 0) {
      return s;
    }
    return reverse;

  }
  public void testMergableDouble() throws IOException {
    final long[][] paths = {{1, 2}, {3, 2}, {2, 4}, {4, 5}, {4, 6}};
    final String[] contigs = {"AAAACC", "CCGGT", "ATCC", "GTAGA", "GACCAC", "GAGTG"};
    final Map <String, String> attr = new HashMap<>();
    attr.put(GraphKmerAttribute.READ_COUNT, "foo");
    attr.put(Consensus.COMBINED, "foo");
    final MutableGraph graph = GraphMapCliTest.makeGraph(2, contigs
        , paths, attr, attr);
    addPath(longs(1, 2, 4, 5), 5, graph);
    addPath(longs(3, 2, 4, 6), 5, graph);
    final MergeNodes merge = new MergeNodes(graph, 2, 2);
    merge.simplifyGraph();
    final StoreDirString dir = new StoreDirString();
    GraphWriter.write(graph, dir, "foo", Collections.<UUID>emptySet());
//    System.err.println(dir);
    final Graph compact = graph.compact();
    final String[] expected = {"AAAACCGGTAGACCAC", "ATCCGGTAGAGTG"};
    finalContigs(expected, compact);
    assertEquals("1/2/4/5", compact.contigAttribute(1, Consensus.COMBINED));
    assertEquals("3/8:(2/4/6)", compact.contigAttribute(2, Consensus.COMBINED));
    assertEquals(0, compact.numberPaths());
  }

  private void newContigs(int oldLength, String[] expected, MutableGraph graph) {
    final Set<String> newContigs = new HashSet<>();
    for (long i = oldLength + 1; i <= graph.numberContigs(); ++i) {
      newContigs.add(ContigString.contigSequenceString(graph.contig(i)));
    }
    final Set<String> expectedContigs = new HashSet<>();
    Collections.addAll(expectedContigs, expected);
    assertEquals(expectedContigs, newContigs);
  }

  private void contigDeleted(MutableGraph graph, List<Long> deleted) {
    for (long i = 1; i <= graph.numberContigs(); ++i) {
      if (deleted.contains(i)) {
        assertTrue("contig " + i + " was not deleted but should have been", graph.contigDeleted(i));
      } else {
        assertFalse("contig " + i + " was deleted but shouldn't have been", graph.contigDeleted(i));
      }
    }
  }

  private void pathDeleted(MutableGraph graph, List<Long> deleted) {
    for (long i = 1; i <= graph.numberPaths(); ++i) {
      if (deleted.contains(i)) {
        assertTrue("path " + i + " was not deleted but should have been", graph.pathDeleted(i));
      } else {
        assertFalse("path " + i + " was deleted but shouldn't have been", graph.pathDeleted(i));
      }
    }
  }

  public void testPalindrome() {
    final long[][] paths = {{1, 2}, {1, -2}, {2, 3}, {-2, 3}, {3, 4}};
    final String[] contigs = {"AAAACC", "CCGG", "GGTT", "TTCCCC", "GTTCCCCC"};
    final Map <String, String> attr = new HashMap<>();
    attr.put(GraphKmerAttribute.READ_COUNT, "foo");
    attr.put(Consensus.COMBINED, "foo");
    final MutableGraph graph = GraphMapCliTest.makeGraph(2, contigs
        , paths, attr, attr);

    MergeNodes.updatePaths(graph, 5, Arrays.asList(3L, 4L));
    assertEquals(2, graph.path(7).contig(0));
    assertEquals(5, graph.path(7).contig(1));
    assertEquals(-2, graph.path(6).contig(0));
    assertEquals(5, graph.path(6).contig(1));

  }
}
