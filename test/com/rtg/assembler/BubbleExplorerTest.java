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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.assembler.DeBruijnGraphBuilderTest.ComparisonNode;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.util.Histogram;

import junit.framework.TestCase;

/**
 */
public class BubbleExplorerTest extends TestCase {
  long[] link(long first, long second) {
    return new long[] {first, second};
  }
  static List<Long> intLinks(long ... links) {
    final List<Long> list = new ArrayList<>();
    for (final long l : links) {
      list.add(l);
    }
    return list;
  }
  static void compareGraph(Map<Long, ComparisonNode> expected, Graph actual) {
    DeBruijnGraphBuilderTest.compareGraph(expected, actual);
  }

  static GraphKmerAttribute buildGraph(int overlap, List<String> contigs, int[] weights, long[]... links) {
    assert contigs.size() == weights.length;
    final GraphKmerAttribute graph = new GraphKmerAttribute(overlap);
    final long[] ids = new long[contigs.size() + 1];
    for (int i = 0; i < contigs.size(); ++i) {
      ids[i + 1] = graph.addContig(new ContigString(contigs.get(i)));
      graph.setKmerFreq(ids[i + 1], weights[i]);
    }
    for (final long[] link : links) {
      assert link.length == 2;
      graph.addPath(new PathArray(convert(link[0], ids), convert(link[1], ids)));
    }
    return graph;
  }
  static long convert(long id, long[] table) {
    final long sign = id > 0 ? 1 : -1;
    return table[(int) Math.abs(id)] * sign;
  }

  GraphKmerAttribute getGraph() {
    return buildGraph(0, Arrays.asList("AA", "AC", "AG", "AT")
        , new int[] {1, 2, 1, 1}
    , link(1, 2)
    , link(1, 3)
    , link(2, 4)
    , link(3, 4)
        );

  }
  public void testSimplestCase() {
    final GraphKmerAttribute graph = getGraph();
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path, 1, 2, 4);

    final BubbleExplorer.Bubble pathReverse = be.explore(-4L);
    checkPath(pathReverse, -4, -2, -1);

    final Histogram histogram = be.popBubbles();
    assertEquals(5, histogram.getLength());
    for (int i = 0; i < 4; ++i) {
      assertEquals(0, histogram.getValue(i));
    }
    assertEquals(1, histogram.getValue(4));
    final Map<Long, ComparisonNode> checker = new HashMap<>();
    checker.put(1L, new ComparisonNode("AAACAT", intLinks(), intLinks()));
    compareGraph(checker, graph);
  }
  public void testStartReversed() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("AA", "AC", "AG", "AT")
        , new int[]{1, 2, 1, 1}
    , link(-1, 2)
    , link(-1, 3)
    , link(2, 4)
    , link(3, 4)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(-1L);
    checkPath(path, -1, 2, 4);

    final BubbleExplorer.Bubble pathReverse = be.explore(-4L);
    checkPath(pathReverse, -4, -2, 1);

    final Histogram histogram = be.popBubbles();
    assertEquals(5, histogram.getLength());
    for (int i = 0; i < 4; ++i) {
      assertEquals(0, histogram.getValue(i));
    }
    assertEquals(1, histogram.getValue(4));
    final Map<Long, ComparisonNode> checker = new HashMap<>();
    checker.put(1L, new ComparisonNode("TTACAT", intLinks(), intLinks()));
    compareGraph(checker, graph);
  }

  void checkPath(BubbleExplorer.Bubble bubble, long ... nodes) {
    if (bubble.mResult != BubbleExplorer.BubbleResult.FOUND) {
      assertEquals(nodes.length, 0);
      return;
    }
    final List<Long> path = bubble.mBestPath;
    final List<Long> nodeList = new ArrayList<>();
    for (final long node : nodes) {
      nodeList.add(node);
    }
    assertEquals(nodeList, path);
  }
  public void testSimpleReverse() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("AA", "AC", "AG", "AT")
        , new int[] {1, 2, 1, 1}
    , link(1, -2)
    , link(1, 3)
    , link(-2, 4)
    , link(3, 4)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    //    System.err.println(path);
    checkPath(path, 1, -2, 4);

    final BubbleExplorer.Bubble path2 = be.explore(-4L);
    checkPath(path2, -4, 2, -1);

    be.popBubbles();
    final Map<Long, DeBruijnGraphBuilderTest.ComparisonNode> checker = new HashMap<>();
    checker.put(1L, new ComparisonNode("AAGTAT", intLinks(), intLinks()));
    compareGraph(checker, graph);
  }
  public void testDeadEnd() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("AA", "AC", "AG", "AT", "AT")
        , new int[] {1, 1, 1, 1, 1}
    , link(1, -2)
    , link(1, 3)
    , link(-4, 2)
    , link(-5, 2)
    , link(3, 4)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path);
    final BubbleExplorer.Bubble path2 = be.explore(-4L);
    checkPath(path2);
  }
  public void testJustUnderMaxLength() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("AA", "ACGTA", "ACGTACGTAC", "AT", "CGTAC")
        , new int[] {1, 2, 2, 1, 1}
    , link(1, -2)
    , link(1, 3)
    , link(-5, 2)
    , link(5, 4)
    , link(3, 4)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path, 1, -2, 5, 4);
  }
  public void testMaxLength() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("AA", "ACGTA", "ACGTACGTAC", "AT", "CGTACG")
        , new int[] {1, 2, 2, 1, 1}
    , link(1, -2)
    , link(1, 3)
    , link(-5, 2)
    , link(5, 4)
    , link(3, 4)
        );

    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path);
  }

  public void testDifference() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("ACA", "AAA", "ACAAAA", "AAA", "AAA", "AAA")
        , new int[] {1, 2, 1, 1, 1, 1}
    , link(1, -2)
    , link(1, 4)
    , link(-2, -3)
    , link(-3, 6)
    , link(4, 5)
    , link(5, 6)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 3, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path);
  }
  public void testFlippedDifference() {
    final GraphKmerAttribute graph = buildGraph(2, Arrays.asList("ACA", "AAA", "AAA", "AAA", "ACAAAA", "AAA")
        , new int[] {1, 2, 1, 1, 1, 1}
    , link(1, -2)
    , link(1, 4)
    , link(-2, -3)
    , link(-3, 6)
    , link(4, 5)
    , link(5, 6)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 3, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path);
  }
  public void testJustUnderDifference() {
    final GraphKmerAttribute graph = buildGraph(2, Arrays.asList("ACA", "AAA", "ACAAA", "AAA", "AAA", "AAA")
        , new int[] {1, 2, 1, 1, 1, 1}
    , link(1, -2)
    , link(1, 4)
    , link(-2, -3)
    , link(-3, 6)
    , link(4, 5)
    , link(5, 6)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 3, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path, 1, -2, -3, 6);
  }
  public void testStartSelfLoop() {
    final GraphKmerAttribute graph = buildGraph(2, Arrays.asList("A", "A", "A", "A")
        , new int[] {1, 2, 1, 1}
    , link(1, 1)
    , link(1, 2)
    , link(1, 3)
    , link(2, 4)
    , link(3, 4)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path);
  }
  public void testStartPreviousLoop() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("A", "A", "A", "A", "A")
        , new int[] {1, 2, 1, 1, 1}
    , link(1, 5)
    , link(1, 2)
    , link(1, 3)
    , link(2, 4)
    , link(3, 4)
    , link(5, 1)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path);
  }
  public void testHairpinEnd() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("A", "A", "A", "A")
        , new int[] {1, 2, 1, 1}
    , link(1, 2)
    , link(1, 3)
    , link(2, 4)
    , link(3, 4)
    , link(4, -4)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path, 1, 2, 4);
  }
  public void testHairpinMiddle() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("A", "A", "A", "A")
        , new int[] {1, 2, 1, 1}
    , link(1, 2)
    , link(1, 3)
    , link(2, 4)
    , link(2, -2)
    , link(3, 4)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path);
  }
  public void testDoubleLoopBack() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("A", "A", "A")
        , new int[] {1, 2, 1}
    , link(1, 2)
    , link(1, 3)
    , link(2, 1)
    , link(3, 1)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path);
  }
  public void testTriplePath() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("A", "C", "G", "A", "C", "G", "G", "T", "G")
        , new int[] {1, 8, 1, 1, 1, 1, 1, 2, 1}
    , link(1, 2)
    , link(2, 3)
    , link(2, 9)
    , link(3, 4)
    , link(3, 8)
    , link(4, 5)
    , link(5, 6)
    , link(6, 7)
    , link(8, 5)
    , link(9, 6)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(2L);
    checkPath(path, 2, 3, 8, 5, 6);

    be.popBubbles();
    //    System.err.println(graph);
    final Map<Long, ComparisonNode> checker = new HashMap<>();
    checker.put(1L, new ComparisonNode("ACGTCGG", intLinks(), intLinks()));
    compareGraph(checker, graph);
  }
  public void testFinalUnexpandable() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("A", "A", "A", "A", "A")
        , new int[] {1, 8, 1, 1, 1}
    , link(1, 2)
    , link(2, 4)
    , link(3, 4)
    , link(5, 4)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    final BubbleExplorer.Bubble path = be.explore(2L);
    checkPath(path);

  }
  public void testBubblesInSeries() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("AA", "AC", "AG", "AT", "CA", "CC", "CG")
        , new int[] {1, 8, 1, 8, 1, 1, 1}
    , link(1, 2)
    , link(1, 6)
    , link(2, 3)
    , link(3, 4)
    , link(3, 7)
    , link(4, 5)
    , link(6, 3)
    , link(7, 5)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    checkPath(be.explore(1L), 1, 2, 3);
    checkPath(be.explore(3L), 3, 4, 5);

    be.popBubbles();
    final Map<Long, ComparisonNode> checker = new HashMap<>();
    checker.put(1L, new ComparisonNode("AAACAGATCA", intLinks(), intLinks()));
    compareGraph(checker, graph);
  }

  public void testNestedBubbles() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("AA", "ACCC", "AG", "AT", "CA", "CC", "CG", "CTTT")
        , new int[] {1, 1, 1, 8, 1, 1, 2, 1}
    , link(1, 2)
    , link(1, 4)
    , link(1, 8)
    , link(2, 3)
    , link(4, 5)
    , link(4, 7)
    , link(5, 6)
    , link(6, 3)
    , link(7, 6)
    , link(8, 3)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    checkPath(be.explore(4L), 4, 7, 6);
    checkPath(be.explore(1L), 1, 4, 7, 6, 3);

    final Histogram histogram = be.popBubbles();
    assertEquals(9, histogram.getLength());
    for (int i = 0; i < 8; ++i) {
      assertEquals(0, histogram.getValue(i));
    }
    assertEquals(1, histogram.getValue(8));
    final Map<Long, ComparisonNode> checker = new HashMap<>();
    checker.put(1L, new ComparisonNode("AAATCGCCAG", intLinks(), intLinks()));
    compareGraph(checker, graph);
  }

  // portion of a failing case from an E.faecalis run
  public void testRealWorldPalindrome() {
    final GraphKmerAttribute graph = buildGraph(0, Arrays.asList("AA", "AC", "AG")
        , new int[] {1, 1, 1}
    , link(1, 3)
    , link(-1, 2)
    , link(2, 1)
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    be.popBubbles();
    //System.err.println(graph);
    final Map<Long, ComparisonNode> checker = new HashMap<>();
    checker.put(1L, new ComparisonNode("AA", intLinks(2, -2), intLinks(3)));
    checker.put(2L, new ComparisonNode("AC", intLinks(-1), intLinks(1)));
    checker.put(3L, new ComparisonNode("AG", intLinks(1), intLinks()));
    compareGraph(checker, graph);
  }
  public void testUnderComplexityThreshold() {
    final long[][] links = new long[200][2];
    final int[] weights = new int[102];
    final List<String> contigs = new ArrayList<>();
    for (int i = 0; i < 102; ++i) {
      weights[i] = 1;
      contigs.add("A");
    }
    weights[6] = 300;
    for (int i = 0; i < 100; ++i) {
      links[i][0] = 1;
      links[i][1] = i + 2;
      links[i + 100][0] = i + 2;
      links[i + 100][1] = 102;
    }

    final GraphKmerAttribute graph = buildGraph(0, contigs
        , weights
        , links
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    checkPath(be.explore(1), 1, 7, 102);
    be.popBubbles();
    //System.err.println(graph);
    final Map<Long, ComparisonNode> checker = new HashMap<>();
    checker.put(1L, new ComparisonNode("AAA", intLinks(), intLinks()));
    compareGraph(checker, graph);

  }
  public void testOverComplexityThreshold() {
    final long[][] links = new long[202][2];
    final int[] weights = new int[103];
    final List<String> contigs = new ArrayList<>();
    for (int i = 0; i < 103; ++i) {
      weights[i] = 1;
      contigs.add("A");
    }
    weights[6] = 300;
    for (int i = 0; i < 101; ++i) {
      links[i][0] = 1;
      links[i][1] = i + 2;
      links[i + 101][0] = i + 2;
      links[i + 101][1] = 103;
    }

    final GraphKmerAttribute graph = buildGraph(0, contigs
        , weights
        , links
        );
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2);
    checkPath(be.explore(1));


  }
  public void testBubbleThreshold() {
    final GraphKmerAttribute graph = getGraph();
    final BubbleExplorer be = new BubbleExplorer(graph, 1, 10, 2, 0.7);
    final BubbleExplorer.Bubble path = be.explore(1L);
    checkPath(path, 1, 2, 4);
    graph.setKmerFreq(2, 69);
    graph.setKmerFreq(3, 30);

    be.popBubbles();
    Map<Long, ComparisonNode> checker = new HashMap<>();
    checker.put(1L, new ComparisonNode("AA", intLinks(), intLinks(2, 3)));
    checker.put(2L, new ComparisonNode("AC", intLinks(1), intLinks(4)));
    checker.put(3L, new ComparisonNode("AG", intLinks(1), intLinks(4)));
    checker.put(4L, new ComparisonNode("AT", intLinks(2, 3), intLinks()));

    graph.setKmerFreq(2, 70);
    graph.setKmerFreq(3, 30);

    be.popBubbles();
    checker = new HashMap<>();
    checker.put(1L, new ComparisonNode("AAACAT", intLinks(), intLinks()));
    compareGraph(checker, graph);
  }

}
