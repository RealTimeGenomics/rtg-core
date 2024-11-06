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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import com.rtg.assembler.graph.PathsIterator;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.mode.DnaUtils;
import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.store.StoreDirString;

import junit.framework.TestCase;

/**
 */
public class ContigCollectorTest extends TestCase {


  GraphKmerAttribute simple(long[]... paths) {
    final GraphKmerAttribute graph = new GraphKmerAttribute(2);
    return simple(graph, paths);
  }
  GraphKmerAttribute simple(GraphKmerAttribute graph, long[]... paths) {
    for (final long[] path : paths) {
      final long first = path[0];
      final long second = path[1];
      final long highest = Math.max(Math.abs(first), Math.abs(second));
      while (graph.numberContigs() < highest) {
        graph.addContig(new ContigString("ACT"));
      }
      graph.addPath(new PathArray(first, second));
    }
    return graph;
  }
  long[] link(long a, long b) {
    return new long[] {a, b};
  }

  public void testUniqueNext() {
    // I don't think we care at this stage about the actual sequences
    final GraphKmerAttribute graph = tipsInMiddleGraph();

    final IntChunks endTips = new IntChunks(graph.numberContigs() + 1);
    final IntChunks startTips = new IntChunks(graph.numberContigs() + 1);
    endTips.set(4, 2);
    startTips.set(5, 2);
    final ContigCollector collector = new ContigCollector(4, 2, startTips, endTips, graph);
    assertEquals(0, collector.uniqueNext(1, false));
    assertEquals(0, collector.uniqueNext(8, true));

    assertEquals(9, collector.uniqueNext(3, true));
    assertEquals(9, collector.uniqueNext(6, false));

    assertEquals(0, collector.uniqueNext(3, false));
    assertEquals(0, collector.uniqueNext(6, true));

    assertEquals(3, collector.uniqueNext(9, false));
    assertEquals(6, collector.uniqueNext(9, true));

    assertEquals(6, collector.uniqueNext(7, false));
    assertEquals(3, collector.uniqueNext(1, true));

  }
  public void testMonogamousTips() {
    // I don't think we care at this stage about the actual sequences
    final GraphKmerAttribute graph = tipsInMiddleGraph();
    final IntChunks endTips = new IntChunks(graph.numberContigs() + 1);
    final IntChunks startTips = new IntChunks(graph.numberContigs() + 1);
    endTips.set(4, 2);
    startTips.set(5, 2);
    final ContigCollector collector = new ContigCollector(4, 2, startTips, endTips, graph);

    // 4 has one node to the right, but 9 has two of which one is non-tip
    assertEquals(9, collector.uniqueNext(4, true));
    assertEquals(0, collector.monogamousLink(4, true));
    assertEquals(3, collector.monogamousLink(9, false));
  }

  /*
  Ascii art for next test!

    4   5 (both 4 & 5 are tips)
     \ /
  2   9   7
   \ / \ /
    3   6
   /     \
  1       8

   */
  private GraphKmerAttribute tipsInMiddleGraph() {
    final GraphKmerAttribute graph = new GraphKmerAttribute(2);
    final long[] linkIds = {
        graph.addContig(new ContigString("ACGTT"))
        ,  graph.addContig(new ContigString("TGGTT"))
        ,  graph.addContig(new ContigString("GTTCCC"))
        ,  graph.addContig(new ContigString("AAACCC"))
        ,  graph.addContig(new ContigString("TGTGGG"))
        ,  graph.addContig(new ContigString("TGTCCGG"))
        ,  graph.addContig(new ContigString("CGGCATG"))
        ,  graph.addContig(new ContigString("CGGCATG"))
        ,  graph.addContig(new ContigString("CCCACGGCTGT"))
    };
    return simple(graph,
        link(linkIds[0], linkIds[2]) , link(linkIds[1], linkIds[2]) , link(linkIds[2], linkIds[8]) , link(linkIds[3], linkIds[8])
        , link(linkIds[8], linkIds[4]) , link(linkIds[8], linkIds[5]) , link(linkIds[5], linkIds[6]) , link(linkIds[5], linkIds[7])
        );
  }


  /*
  Ascii art for next test!

  3 (1, 2 & 4 are tips)
   \
  2-4--5
   /
  1

   */
  public void testSingleTip() {
    // I don't think we care at this stage about the actual sequences
    final GraphKmerAttribute graph = simple(
        link(1, 4) , link(2, 4) , link(3, 4), link(4, 5)
        );

    final IntChunks endTips = new IntChunks(graph.numberContigs() + 1);
    final IntChunks startTips = new IntChunks(graph.numberContigs() + 1);
    endTips.set(1, 2);
    startTips.set(2, 2);
    startTips.set(3, 2);
    startTips.set(5, 2);
    final ContigCollector collector = new ContigCollector(4, 2, startTips, endTips, graph);
    assertEquals(0, collector.uniqueNext(4, false));
    assertEquals(5, collector.uniqueNext(4, true));

  }

  public void testRcNext() {
    // I don't think we care at this stage about the actual sequences
    final GraphKmerAttribute graph = reversedNodeGraph();
    final IntChunks endTips = new IntChunks(graph.numberContigs() + 1);
    final IntChunks startTips = new IntChunks(graph.numberContigs() + 1);
    final ContigCollector collector = new ContigCollector(4, 2, startTips, endTips, graph);
    assertEquals(-4, collector.uniqueNext(3, true));
    assertEquals(-4, collector.uniqueNext(5, false));

    assertEquals(-3, collector.uniqueNext(4, true));
    assertEquals(-5, collector.uniqueNext(4, false));

    assertEquals(3, collector.uniqueNext(-4, false));
    assertEquals(5, collector.uniqueNext(-4, true));
    assertEquals(0, collector.uniqueNext(5, true));
    assertEquals(0, collector.uniqueNext(-5, false));

  }

  /*
  Ascii art for next graph!  No tips here

  2            6
   \          /
    3--(-4)--5
   /          \
  1            (-7)

   */
  private GraphKmerAttribute reversedNodeGraph() {
    return simple(
        link(1, 3) , link(2, 3) , link(3, -4) , link(-5, 4)
        , link(5, 6) , link(5, -7)
        );
  }

  public void testMonogamous() {
    final GraphKmerAttribute graph = reversedNodeGraph();
    final IntChunks endTips = new IntChunks(graph.numberContigs() + 1);
    final IntChunks startTips = new IntChunks(graph.numberContigs() + 1);
    final ContigCollector collector = new ContigCollector(4, 2, startTips, endTips, graph);
    assertEquals(0, collector.monogamousLink(5, true));
    assertEquals(0, collector.monogamousLink(-7, true));
  }

  public void testWalk() {
    final GraphKmerAttribute graph = reversedNodeGraph();
    final IntChunks endTips = new IntChunks(graph.numberContigs() + 1);
    final IntChunks startTips = new IntChunks(graph.numberContigs() + 1);
    final ContigCollector collector = new ContigCollector(4, 2, startTips, endTips, graph);
    final List<Long> walked = new ArrayList<>();

    walked.add(3L);
    collector.walk(walked, 3, true);
    assertEquals(Arrays.asList(3L, -4L, 5L), walked);

    walked.clear();
    walked.add(3L);
    collector.walk(walked, 3, false);
    assertEquals(Arrays.asList(3L), walked);

    walked.clear();
    walked.add(-4L);
    collector.walk(walked, -4, true);
    assertEquals(Arrays.asList(-4L, 5L), walked);
    collector.walk(walked, -4, false);
    assertEquals(Arrays.asList(3L, -4L, 5L), walked);
  }

  public void testCollect() {
    final GraphKmerAttribute graph = tipsInMiddleGraph();
    final IntChunks endTips = new IntChunks(graph.numberContigs() + 1);
    final IntChunks startTips = new IntChunks(graph.numberContigs() + 1);
    endTips.set(4, 2);
    startTips.set(5, 2);
    graph.setKmerFreq(3, 3);
    graph.setKmerFreq(9, 2);
    graph.setKmerFreq(6, 1);
    final ContigCollector collector = new ContigCollector(4, 4, startTips, endTips, graph);
    collector.collapse();
    assertTrue(graph.contigDeleted(3));
    assertTrue(graph.contigDeleted(9));
    assertTrue(graph.contigDeleted(6));
    final int combinedId = 10;
    assertEquals(combinedId, graph.numberContigs());
    final String expectedContig = "GTTCCCACGGCTGTCCGG";
    final String actualContig = ContigString.contigSequenceString(graph.contig(combinedId));
    assertTrue("expected <" + expectedContig + "> or it's RC but was <" + actualContig + ">", expectedContig.equals(actualContig) || DnaUtils.reverseComplement(expectedContig).equals(actualContig));
    assertTrue(hasPath(graph, combinedId, 1, false));
    assertTrue(hasPath(graph, combinedId, 2, false));
    assertTrue(hasPath(graph, combinedId, 7, true));
    assertTrue(hasPath(graph, combinedId, 8, true));
    assertEquals(4, numberPaths(graph, combinedId));
    assertEquals(6, graph.kmerFreq(combinedId));
  }
  public void testContrivedPastTipArray() {

    final GraphKmerAttribute graph = new GraphKmerAttribute(2);
    graph.addContig(new ContigString("CCCACG"));
    graph.addContig(new ContigString("ACGCCC"));
    graph.addPath(new PathArray(1, 2));
    final IntChunks endTips = new IntChunks(graph.numberContigs());
    final IntChunks startTips = new IntChunks(graph.numberContigs());
    endTips.set(1, 2);
    startTips.set(1, 2);
    final ContigCollector collector = new ContigCollector(4, 4, startTips, endTips, graph);
    collector.collapse();
    final int combinedId = 3;
    assertEquals(combinedId, graph.numberContigs());
    final String expectedContig = "CCCACGCCC";
    final String actualContig = ContigString.contigSequenceString(graph.contig(combinedId));
    assertTrue("expected <" + expectedContig + "> or it's RC but was <" + actualContig + ">", expectedContig.equals(actualContig) || DnaUtils.reverseComplement(expectedContig).equals(actualContig));
    assertEquals(0, numberPaths(graph, combinedId));
    assertEquals(0, startTips.getInt(combinedId));
    assertEquals(0, endTips.getInt(combinedId));

  }
  public void testTipMerging() {

    final GraphKmerAttribute graph = new GraphKmerAttribute(2);
    graph.addContig(new ContigString("CCCACG"));
    graph.addContig(new ContigString("ACGCCC"));
    graph.addPath(new PathArray(1, 2));
    final IntChunks endTips = new IntChunks(graph.numberContigs() + 1);
    final IntChunks startTips = new IntChunks(graph.numberContigs() + 1);
    endTips.set(1, 2);
    endTips.set(2, 5);
    final ContigCollector collector = new ContigCollector(4, 4, startTips, endTips, graph);
    collector.collapse();
    final int combinedId = 3;
    assertEquals(combinedId, graph.numberContigs());
    final String expectedContig = "CCCACGCCC";
    final String actualContig = ContigString.contigSequenceString(graph.contig(combinedId));
    assertTrue("expected <" + expectedContig + "> or it's RC but was <" + actualContig + ">", expectedContig.equals(actualContig) || DnaUtils.reverseComplement(expectedContig).equals(actualContig));
    assertEquals(0, numberPaths(graph, combinedId));
    assertEquals(endTips.getInt(combinedId), 5);
  }

  public void testTipMergingReverse() {

    final GraphKmerAttribute graph = new GraphKmerAttribute(2);
    graph.addContig(new ContigString("CGTGGG"));
    graph.addContig(new ContigString("ACGCCC"));
    graph.addPath(new PathArray(-1, 2));
    final IntChunks endTips = new IntChunks(graph.numberContigs() + 1);
    final IntChunks startTips = new IntChunks(graph.numberContigs() + 1);
    endTips.set(1, 5);
    endTips.set(2, 7);
    startTips.set(2, 9);
    final ContigCollector collector = new ContigCollector(4, 4, startTips, endTips, graph);
    collector.collapse();
    final int combinedId = 3;
    assertEquals(combinedId, graph.numberContigs());
    final String expectedContig = "CCCACGCCC";
    final String actualContig = ContigString.contigSequenceString(graph.contig(combinedId));
    assertTrue("expected <" + expectedContig + "> or it's RC but was <" + actualContig + ">", expectedContig.equals(actualContig) || DnaUtils.reverseComplement(expectedContig).equals(actualContig));
    assertEquals(0, numberPaths(graph, combinedId));

    assertEquals(0, startTips.getInt(combinedId));
    assertEquals(9, endTips.getInt(combinedId));
  }
  public void testTipMergingReverse2() {

    final GraphKmerAttribute graph = new GraphKmerAttribute(2);
    graph.addContig(new ContigString("CGTGGG"));
    graph.addContig(new ContigString("ACGCCC"));
    graph.addPath(new PathArray(-1, 2));
    final IntChunks endTips = new IntChunks(graph.numberContigs() + 1);
    final IntChunks startTips = new IntChunks(graph.numberContigs() + 1);
    endTips.set(1, 2);
    endTips.set(2, 7);
    startTips.set(2, 9);
    final ContigCollector collector = new ContigCollector(4, 4, startTips, endTips, graph);
    collector.collapse();
    final int combinedId = 3;
    assertEquals(combinedId, graph.numberContigs());
    final String expectedContig = "CCCACGCCC";
    final String actualContig = ContigString.contigSequenceString(graph.contig(combinedId));
    assertTrue("expected <" + expectedContig + "> or it's RC but was <" + actualContig + ">", expectedContig.equals(actualContig) || DnaUtils.reverseComplement(expectedContig).equals(actualContig));
    assertEquals(0, numberPaths(graph, combinedId));

    assertEquals(0, endTips.getInt(combinedId));
    assertEquals(9, startTips.getInt(combinedId));
  }
  int numberPaths(GraphKmerAttribute graph, long id) {
    final PathsIterator pathsIterator = graph.paths(id);
    int count = 0;
    while (pathsIterator.nextPathId() != 0) {
      ++count;
    }
    return count;
  }
  boolean hasPath(GraphKmerAttribute graph, long id, long linkId, boolean end) {
    final PathsIterator pathsIterator = graph.paths(id);
    long pathId;
    boolean found = false;
    while ((pathId = pathsIterator.nextPathId()) != 0) {
      final int index = pathsIterator.contigIndex();
      final long linkedContig = graph.pathContig(pathId, 1 - index);
      if (end && index == 0 || !end && index == 1) {
        if (linkedContig == linkId) {
          found = true;
          break;
        }
      }
    }
    return found;
  }

  public void testCollapseWithLongPaths() throws IOException {
    final HashMap<String, String> pathAttr = new HashMap<>();
    pathAttr.put("readCount", "count em all");
    final GraphKmerAttribute g = GraphMapCliTest.makeGraph(0, new String[]{"AAAA", "GGGG", "AAA", "A", "AA", "GGGGG"}, new long[][]{{1, 2, 3, 4}, {5, 2, 3, 6}}, Collections.emptyMap(), pathAttr);
    g.setPathAttribute(1, "readCount", "10");
    g.setPathAttribute(2, "readCount", "10");
    new ContigCollector(0, 1, null, null, g).collapse();
    final StoreDirString str = new StoreDirString();
    GraphWriter.writeWithDeleted(g, str, "foo", Collections.emptySet());
    assertEquals("GGGGAAA", ContigString.contigSequenceString(g.contig(7)));
    assertEquals(Arrays.asList(1L, 7L, 4L), MergeNodes.getPath(g, 4));
    assertEquals(Arrays.asList(5L, 7L, 6L), MergeNodes.getPath(g, 3));
  }

  public void testCollapseWithLongerPaths() throws IOException {
    final HashMap<String, String> pathAttr = new HashMap<>();
    pathAttr.put("readCount", "count em all");
    final HashMap<String, String> contigAttr = new HashMap<>();
    contigAttr.put(Consensus.COMBINED, "count em all");
    final GraphKmerAttribute g = GraphMapCliTest.makeGraph(0,
        new String[]{"AAAA", "GGGG", "AAA", "A", "AA", "GGGGG", "T", "C"}
        , new long[][]{{1, 2, 3, 4}, {5, 2, 3, 6}, {7, 1, 2, 3, 4}, {8, 1, 2, 3, 4}}
        , contigAttr, pathAttr);
    g.setPathAttribute(1, "readCount", "10");
    g.setPathAttribute(2, "readCount", "10");
    new ContigCollector(0, 1, null, null, g).collapse();
    final StoreDirString str = new StoreDirString();
    GraphWriter.writeWithDeleted(g, str, "foo", Collections.emptySet());
    assertEquals("GGGGAAA", ContigString.contigSequenceString(g.contig(9)));
    assertEquals(Arrays.asList(1L, 9L, 4L), MergeNodes.getPath(g, 8));
    assertEquals(Arrays.asList(5L, 9L, 6L), MergeNodes.getPath(g, 7));
    assertEquals(Arrays.asList(7L, 1L, 9L, 4L), MergeNodes.getPath(g, 6));
    assertEquals(Arrays.asList(8L, 1L, 9L, 4L), MergeNodes.getPath(g, 5));
    assertEquals("2/3", g.contigAttribute(9, Consensus.COMBINED));

    assertEquals("7/1/-9:(2/3)", ContigCollector.buildCombinedString(Arrays.asList(7L, 1L, -9L), g));
  }
}
