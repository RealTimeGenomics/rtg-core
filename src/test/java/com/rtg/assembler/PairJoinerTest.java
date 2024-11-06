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
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;

import junit.framework.TestCase;

/**
 */
public class PairJoinerTest extends TestCase {
  public void testSingleContig() {
    final Graph graph = GraphMapCliTest.makeGraph(0, new String[]{"ACGTTTTTTTGGTTG"}, new long[][]{});
    final PairJoiner joiner = new PairJoiner(graph, 10);
    Set<GraphAlignment> destinations = new HashSet<>();
    destinations.add(new GraphAlignment(0, 4, Arrays.asList(-1L), 2, graph));
    final GraphAlignment start = new GraphAlignment(0, 3, Arrays.asList(1L), 1, graph);
    // FR
    Set<GraphAlignment> actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 6, 7);
    List<GraphAlignment> expected = Arrays.asList(new GraphAlignment(0, 14, Arrays.asList(1L), 3, graph));
    checkAlignments(expected, actual);

    // Too far apart
    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 7, 8);
    expected = Collections.emptyList();
    checkAlignments(expected, actual);

    // Too close together
    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 5, 6);
    expected = Collections.emptyList();
    checkAlignments(expected, actual);

    // FF
    destinations = new HashSet<>();
    destinations.add(new GraphAlignment(10, 14, Arrays.asList(1L), 2, graph));
    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 5, 6);
    expected = Collections.emptyList();
    checkAlignments(expected, actual);

    // RF
    actual = joiner.joinAlignments(new GraphAlignment(11, 14, Arrays.asList(-1L), 1, graph), PairJoiner.destinationMap(destinations, graph), 5, 6);
    expected = Collections.emptyList();
    checkAlignments(expected, actual);
  }

  public void testX() {
    // 12345678901234
    //  --123456----
    // ACGT    TGGTTG
    //    TTTTTT
    // CCCT    TCCACCG
    //    12345-------
    final Graph graph = getXGraph();
    final PairJoiner joiner = new PairJoiner(graph, 1);

    Set<GraphAlignment> destinations = new HashSet<>();
    final GraphAlignment destination1 = new GraphAlignment(1, 4, Arrays.asList(-3L), 2, graph);
    destinations.add(destination1);
    final GraphAlignment destination2 = new GraphAlignment(0, 6, Arrays.asList(-5L), 4, graph);
    destinations.add(destination2);

    final GraphAlignment start = new GraphAlignment(1, 2, Arrays.asList(1L), 1, graph);
    Set<GraphAlignment> actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 5, 7);
    final GraphAlignment alignment1 = new GraphAlignment(1, 4, Arrays.asList(1L, 2L, 3L), 3, graph);
    final GraphAlignment alignment2 = new GraphAlignment(1, 6, Arrays.asList(1L, 2L, 5L), 5, graph);
    List<GraphAlignment> expected = Arrays.asList(alignment1
    , alignment2);
    checkAlignments(expected, actual);

    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 6, 8);
    expected = Arrays.asList(alignment1);
    checkAlignments(expected, actual);

    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 4, 6);
    expected = Arrays.asList(alignment2);
    checkAlignments(expected, actual);

    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 3, 5);
    expected = Collections.emptyList();
    checkAlignments(expected, actual);

    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 7, 9);
    expected = Collections.emptyList();
    checkAlignments(expected, actual);


    destinations = new HashSet<>();
    destinations.add(start);
    destinations.add(destination2);

    actual = joiner.joinAlignments(destination1, PairJoiner.destinationMap(destinations, graph), 5, 7);
    expected = Arrays.asList(new GraphAlignment(1, 2, Arrays.asList(-3L, -2L, -1L), 3, graph));
    checkAlignments(expected, actual);

    destinations = new HashSet<>();
    destinations.add(start);
    destinations.add(destination1);

    actual = joiner.joinAlignments(destination2, PairJoiner.destinationMap(destinations, graph), 5, 7);
    expected = Arrays.asList(new GraphAlignment(0, 2, Arrays.asList(-5L, -2L, -1L), 5, graph));
    checkAlignments(expected, actual);
  }

  private void checkAlignments(List<GraphAlignment> expected, Collection<GraphAlignment> actual) {
    assertTrue("Expected <" + expected + "> but was <" + actual + ">", actual.size() == expected.size() && actual.containsAll(expected));
  }

  public void testLoop() {

    final Graph graph = GraphMapCliTest.makeGraph(0, new String[]{"ACGT", "TTTTTT", "TGGTTG"}, new long[][]{{1, 2}, {2, 3}, {2, 2}});
    final PairJoiner joiner = new PairJoiner(graph, 1);
    final Set<GraphAlignment> destinations = new HashSet<>();
    final GraphAlignment destination = new GraphAlignment(1, 4, Arrays.asList(-3L), 2, graph);
    destinations.add(destination);
    GraphAlignment start = new GraphAlignment(1, 2, Arrays.asList(1L), 1, graph);
    final GraphAlignment oneLoop = new GraphAlignment(1, 4, Arrays.asList(1L, 2L, 3L), 3, graph);
    Set<GraphAlignment> actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 6, 7);
    List<GraphAlignment> expected = Arrays.asList(oneLoop);
    checkAlignments(expected, actual);

    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 11, 12);
    final GraphAlignment twoLoops = new GraphAlignment(1, 4, Arrays.asList(1L, 2L, 2L, 3L), 3, graph);
    expected = Arrays.asList(twoLoops);
    checkAlignments(expected, actual);

    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 6, 12);
    expected = Arrays.asList(oneLoop, twoLoops);
    checkAlignments(expected, actual);

    // Aligning RF on a looping contig
    destinations.clear();
    destinations.add(new GraphAlignment(3, 4, Arrays.asList(-2L), 2, graph));
    start = new GraphAlignment(2, 3, Arrays.asList(2L), 1, graph);
    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 2, 3);
    expected = Arrays.asList(new GraphAlignment(2, 2, Arrays.asList(2L, 2L), 3, graph));
    checkAlignments(expected, actual);

    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 7, 8);
    expected = Arrays.asList(new GraphAlignment(2, 2, Arrays.asList(2L, 2L, 2L), 3, graph));
    checkAlignments(expected, actual);
  }
  public void testGodDamnedPalindrome() {
    final Graph graph = GraphMapCliTest.makeGraph(0, new String[]{"ACGT", "TTTTTT", "TATA"}, new long[][]{{1, 2}, {2, 3}, {2, -3}});
    final PairJoiner joiner = new PairJoiner(graph, 1);

    final GraphAlignment destination = new GraphAlignment(1, 4, Arrays.asList(2L), 2, graph);
    final GraphAlignment start = new GraphAlignment(1, 2, Arrays.asList(1L), 1, graph);
    final Set<GraphAlignment> destinations = new HashSet<>();
    destinations.add(destination);

    final Set<GraphAlignment> actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 9, 10);
    final List<GraphAlignment> expected = Arrays.asList(
        new GraphAlignment(1, 4, Arrays.asList(1L, 2L, 3L, -2L), 3, graph)
    );
    checkAlignments(expected, actual);
  }

  public void testGodDamnedPalindromeEnd() {
    final Graph graph = GraphMapCliTest.makeGraph(0, new String[]{"ACCT", "TTTTTT", "TATA", "GGGGG"}, new long[][]{{1, 2}, {2, 3}, {2, -3}, {3, 4}, {-3, 4}});
    final PairJoiner joiner = new PairJoiner(graph, 1);

    final GraphAlignment destination = new GraphAlignment(1, 2, Arrays.asList(-4L, -3L), 2, graph);
    final GraphAlignment start = new GraphAlignment(1, 2, Arrays.asList(1L), 1, graph);
    Set<GraphAlignment> destinations = new HashSet<>();
    destinations.add(destination);

    Set<GraphAlignment> actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 5, 20);
    List<GraphAlignment> expected = Arrays.asList(
        new GraphAlignment(1, 3, Arrays.asList(1L, 2L, 3L, 4L), 3, graph)
    );
    checkAlignments(expected, actual);

    final GraphAlignment destination2 = new GraphAlignment(1, 2, Arrays.asList(-4L, 3L), 2, graph);
    destinations = new HashSet<>();
    destinations.add(destination2);
    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 5, 20);
    expected = Arrays.asList(
        new GraphAlignment(1, 3, Arrays.asList(1L, 2L, -3L, 4L), 3, graph)
    );
    checkAlignments(expected, actual);
  }

  public void testComplexEnds() {
    final Graph graph = GraphMapCliTest.makeGraph(0, new String[]{"ACGT", "TTTT", "TATA", "ACCA", "AGGT", "TCCCGGGT", "GAGA"}, new long[][]{{1, 2}, {2, 3}, {3, 4}, {4, 5}, {1, 6}, {7, 5}, {6, 7}});
    final PairJoiner joiner = new PairJoiner(graph, 1);

    final GraphAlignment start = new GraphAlignment(1, 2, Arrays.asList(1L, 2L), 1, graph);
    final GraphAlignment destination = new GraphAlignment(1, 3, Arrays.asList(-5L, -4L), 2, graph);
    final Set<GraphAlignment> destinations = new HashSet<>();
    destinations.add(destination);

    final Set<GraphAlignment> actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), 3, 40);
    final List<GraphAlignment> expected = Arrays.asList(
        new GraphAlignment(1, 2, Arrays.asList(1L, 2L, 3L, 4L, 5L), 3, graph)
    );
    checkAlignments(expected, actual);
  }
  public void testWrapperMethod() {
    final Graph graph = getXGraph();
    final PairJoiner joiner = new PairJoiner(graph, 1);
    final Set<GraphAlignment> starts = new HashSet<>();
    final Set<GraphAlignment> destinations = new HashSet<>();
    starts.add(new GraphAlignment(1, 2, Arrays.asList(1L), 1, graph));
    starts.add(new GraphAlignment(1, 2, Arrays.asList(4L), 1, graph));
    destinations.add(new GraphAlignment(1, 4, Arrays.asList(-3L), 2, graph));
    destinations.add(new GraphAlignment(0, 6, Arrays.asList(-5L), 4, graph));
    final ArrayList<Set<GraphAlignment>> sets = new ArrayList<>();
    sets.add(starts);
    sets.add(destinations);
    final Set<GraphAlignment> actual = joiner.paired(sets, 0, 100);
    final List<GraphAlignment> expected = Arrays.asList(
        new GraphAlignment(1, 4, Arrays.asList(1L, 2L, 3L), 3, graph)
        , new GraphAlignment(1, 4, Arrays.asList(4L, 2L, 3L), 3, graph)
        , new GraphAlignment(1, 6, Arrays.asList(1L, 2L, 5L), 5, graph)
        , new GraphAlignment(1, 6, Arrays.asList(4L, 2L, 5L), 5, graph)
    );
    checkAlignments(expected, actual);


  }

  private MutableGraph getXGraph() {
    return GraphMapCliTest.makeGraph(0, new String[]{"ACGT", "TTTTTT", "TGGTTG", "CCCT", "TCCACCG"}, new long[][]{{1, 2}, {2, 3}, {4, 2}, {2, 5}});
  }

  public void testHashCodeOfPairChain() {
    // This is necessary to prevent non determinism because we put these in a set
    final PairJoiner.PairChain pc1 = new PairJoiner.PairChain(null, 100, 100);
    final PairJoiner.PairChain pc2 = new PairJoiner.PairChain(null, 100, 100);
    assertEquals(pc1.hashCode(), pc2.hashCode());
    assertTrue(pc1.equals(pc2));
    assertTrue(pc2.equals(pc1));
    assertTrue(pc2.equals(pc2));
  }

  public void testSlopTigsMultiOverlap() {
    final Graph graph = GraphMapCliTest.makeGraph(0,
        new String[]{"ACGTTGTGTTAGAGTGATG"
            , "TAGAGTGATGAAAC"
            , "GTGATGAAACGTTACCCA"
            // -------------------|
            //    |----------------------
        }
        , new long[][]{{1, 2}, {2, 3}}
    );
    final PairJoiner joiner = new PairJoiner(graph, 10);
    final Set<GraphAlignment> starts = new HashSet<>();
    final Set<GraphAlignment> destinations = new HashSet<>();
    starts.add(new GraphAlignment(1, 10, Arrays.asList(1L, 2L, 3L), 1, graph));
    destinations.add(new GraphAlignment(1, 11, Arrays.asList(-3L, -2L, -1L), 1, graph));
    final ArrayList<Set<GraphAlignment>> sets = new ArrayList<>();
    sets.add(starts);
    sets.add(destinations);
    Set<GraphAlignment> actual = joiner.paired(sets, -16, 100);
    final List<GraphAlignment> expected = Arrays.asList(
        new GraphAlignment(1, 16, Arrays.asList(1L, 2L, 3L), 2, graph)
    );
    checkAlignments(expected, actual);
    actual = joiner.paired(sets, -15, 100);
    checkAlignments(Collections.emptyList(), actual);
  }
  public void testSlopTigsKmerOverlap() {
    final Graph graph = GraphMapCliTest.makeGraph(0,
        new String[]{"ACGTTGTGTTAGAGTGATG"
            , "TAGAGTGATGAAAC"
            , "GTGATGAAACGTTACCCA"
            // -------------------|
            //                       |----------------------
        }
        , new long[][]{{1, 2}, {2, 3}}
    );
    final PairJoiner joiner = new PairJoiner(graph, 10);
    final Set<GraphAlignment> starts = new HashSet<>();
    final Set<GraphAlignment> destinations = new HashSet<>();
    starts.add(new GraphAlignment(1, 17, Arrays.asList(1L), 1, graph));
    destinations.add(new GraphAlignment(1, 16, Arrays.asList(-3L), 1, graph));
    final ArrayList<Set<GraphAlignment>> sets = new ArrayList<>();
    sets.add(starts);
    sets.add(destinations);
    Set<GraphAlignment> actual = joiner.paired(sets, -4, 100);
    final List<GraphAlignment> expected = Arrays.asList(
        new GraphAlignment(1, 16, Arrays.asList(1L, 2L, 3L), 2, graph)
    );
    checkAlignments(expected, actual);
    actual = joiner.paired(sets, -3, 100);
    checkAlignments(Collections.emptyList(), actual);
  }

  public void testSlopTigsSameContig() {
    final Graph graph = GraphMapCliTest.makeGraph(0,
        // ---------------|
        //    |---------------
        new String[]{"ACGTTGTGTTAGAGTGATG"
        }
        , new long[][]{}
    );
    final PairJoiner joiner = new PairJoiner(graph, 10);
    final Set<GraphAlignment> starts = new HashSet<>();
    final Set<GraphAlignment> destinations = new HashSet<>();
    starts.add(new GraphAlignment(1, 15, Arrays.asList(1L), 1, graph));
    destinations.add(new GraphAlignment(1, 15, Arrays.asList(-1L), 1, graph));
    final ArrayList<Set<GraphAlignment>> sets = new ArrayList<>();
    sets.add(starts);
    sets.add(destinations);
    Set<GraphAlignment> actual = joiner.paired(sets, -13, 100);
    final List<GraphAlignment> expected = Arrays.asList(
        new GraphAlignment(1, 17, Arrays.asList(1L), 2, graph)
    );
    checkAlignments(expected, actual);
    actual = joiner.paired(sets, -12, 100);
    checkAlignments(Collections.emptyList(), actual);
  }
  public void testSlopTigsAltBranch() {
    final Graph graph = GraphMapCliTest.makeGraph(0,
        new String[]{"ACGTTGTGTTAGAGTGATG"
            , "TAGAGTGATGAAAC"
            , "TAGAGTGATGCCGA"
            // -------------------|
            //                       |----------------------
        }
        , new long[][]{{1, 2}, {1, 3}}
    );
    final PairJoiner joiner = new PairJoiner(graph, 10);
    final Set<GraphAlignment> starts = new HashSet<>();
    final Set<GraphAlignment> destinations = new HashSet<>();
    starts.add(new GraphAlignment(1, 12, Arrays.asList(1L, 2L), 1, graph));
    destinations.add(new GraphAlignment(1, 12, Arrays.asList(-3L, -1L), 1, graph));
    final ArrayList<Set<GraphAlignment>> sets = new ArrayList<>();
    sets.add(starts);
    sets.add(destinations);
    final Set<GraphAlignment> actual = joiner.paired(sets, -100, 100);
    final List<GraphAlignment> expected = Collections.emptyList();
    checkAlignments(expected, actual);
  }
  public void testSloptigPalindrome() {
    final Graph graph = GraphMapCliTest.makeGraph(0,
        new String[]{"ACCCC", "CCCGAC", "ACGT", "GTTTTT", "TTCC"
        }
        , new long[][]{{1, 2}, {2, 3}, {2, -3}, {3, 4}, {-3, 4}, {4, 5}}
    );
    final PairJoiner joiner = new PairJoiner(graph, 1);

    final GraphAlignment destination = new GraphAlignment(0, 2, Arrays.asList(-5L, -4L, 3L, -2L, -1L), 1, graph);
    final GraphAlignment start = new GraphAlignment(1, 3, Arrays.asList(1L, 2L, 3L, 4L, 5L), 1, graph);
    Set<GraphAlignment> destinations = new HashSet<>();
    destinations.add(destination);

    Set<GraphAlignment> actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), -30, 20);
    List<GraphAlignment> expected = Arrays.asList(
        new GraphAlignment(1, 3, Arrays.asList(1L, 2L, -3L, 4L, 5L), 2, graph)
    );
    checkAlignments(expected, actual);

    final GraphAlignment destination2 = new GraphAlignment(0, 2, Arrays.asList(-5L, -4L, -3L, -2L, -1L), 1, graph);
    destinations = new HashSet<>();
    destinations.add(destination2);
    actual = joiner.joinAlignments(start, PairJoiner.destinationMap(destinations, graph), -30, 20);
    expected = Arrays.asList(
        new GraphAlignment(1, 3, Arrays.asList(1L, 2L, 3L, 4L, 5L), 2, graph)
    );
    checkAlignments(expected, actual);
  }
}
