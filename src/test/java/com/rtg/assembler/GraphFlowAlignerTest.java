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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.rtg.assembler.graph.Graph;
import com.rtg.mode.DnaUtils;
import com.rtg.util.IntegerOrPercentage;

import junit.framework.TestCase;

/**
 */
public class GraphFlowAlignerTest extends TestCase {
  Set<GraphAlignment> set(GraphAlignment... alignments) {
    final Set<GraphAlignment> expected = new HashSet<>();
    Collections.addAll(expected, alignments);
    return expected;
  }
  List<Long> longs(long... l) {
    final List<Long> expected = new ArrayList<>();
    for (long g : l) {
      expected.add(g);
    }
    return expected;
  }

  public void testFlowAligner() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"ACGTAAAACTGAAACCCTAAACC"}, new long[][]{});
    final GraphFlowAligner aligner = new GraphFlowAligner(graph, new IntegerOrPercentage(10), new GraphTraversions(graph));
    final Set<GraphAlignment> alignments = aligner.align(DnaUtils.encodeString("AAACTGGGGAACCTAAA"), 0, new ContigPosition(1, 4, graph));
    final Set<GraphAlignment> expected = set(new GraphAlignment(4, 20, longs(1), 6, graph));
    assertEquals(expected, alignments);
  }

  public void testMismatchAligner() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"ACGTAAAACTGAAACCCTAAACC"}, new long[][]{});
    final GraphFlowAligner aligner = new GraphFlowAligner(graph, new IntegerOrPercentage(10), new GraphTraversions(graph));
    final Set<GraphAlignment> alignments = aligner.align(DnaUtils.encodeString("AAAACTCAAACCCTAAACC"), 0, new ContigPosition(1, 4, graph));
    final Set<GraphAlignment> expected = set(new GraphAlignment(4, 22, longs(1), 1, graph));
    assertEquals(expected, alignments);
  }

  public void testAlignmentScore() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"ACGTAAAACTGAAACCCTAAACC"}, new long[][]{});
    GraphFlowAligner aligner = new GraphFlowAligner(graph, new IntegerOrPercentage(4), new GraphTraversions(graph));
    Set<GraphAlignment> alignments = aligner.align(DnaUtils.encodeString("AAAACTGTTTCCCTAAACC"), 0, new ContigPosition(1, 4, graph));
    Set<GraphAlignment> expected = set();
    assertEquals(expected, alignments);

    aligner = new GraphFlowAligner(graph, new IntegerOrPercentage(5), new GraphTraversions(graph));
    alignments = aligner.align(DnaUtils.encodeString("AAAACTGTTTCCCTAAACC"), 0, new ContigPosition(1, 4, graph));
    expected = set(new GraphAlignment(4, 22, longs(1), 5, graph));
    assertEquals(expected, alignments);

  }
  public void testScoreCrossContig() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"ACGTAAAA", "AACTGAACCCTAAAG", "AACTGGGGACCCCTTAAG"}, new long[][]{{1, 2}, {1, 3}});
    final GraphFlowAligner aligner = new GraphFlowAligner(graph, new IntegerOrPercentage(5), new GraphTraversions(graph));
    final Set<GraphAlignment> alignments = aligner.align(DnaUtils.encodeString("AAACTGGGGAACCTAAA"), 0, new ContigPosition(1, 4, graph));
    final Set<GraphAlignment> expected = set(
        new GraphAlignment(4, 13, longs(1, 2), 5, graph)
    );
    assertEquals(expected, alignments);
  }

  public void testFlowAlignerCrossContig() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"ACGTAAAA", "AACTGAACCCTAAAG", "AACTGGGGACCCCTTAAG"}, new long[][]{{1, 2}, {1, 3}});
    final GraphFlowAligner aligner = new GraphFlowAligner(graph, new IntegerOrPercentage(15), new GraphTraversions(graph));
    Set<GraphAlignment> alignments = aligner.align(DnaUtils.encodeString("AAACTGGGGAACCTAAA"), 0, new ContigPosition(1, 4, graph));
    Set<GraphAlignment> expected = set(
        new GraphAlignment(4, 13, longs(1, 2), 5, graph)
        , new GraphAlignment(4, 16, longs(1, 3), 6, graph)
    );
    assertEquals(expected, alignments);

    alignments = aligner.align(DnaUtils.encodeString("AAACTGGGGAACCTAAA"), 2, new ContigPosition(1, 4, graph)); expected = set(
        new GraphAlignment(4, 13, longs(1, 2), 9, graph)
        , new GraphAlignment(4, 16, longs(1, 3), 10, graph)
    );
    assertEquals(expected, alignments);

    alignments = aligner.align(DnaUtils.encodeString("AAACTGGGGAACCTAAA"), 5, new ContigPosition(3, 4, graph)); expected = set(
         new GraphAlignment(5, 16, longs(1, 3), 5, graph)
    );
    assertEquals(expected, alignments);
  }
}
