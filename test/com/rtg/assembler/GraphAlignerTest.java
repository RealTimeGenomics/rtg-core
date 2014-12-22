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
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.rtg.assembler.graph.Graph;
import com.rtg.mode.DnaUtils;
import com.rtg.util.IntegerOrPercentage;

import junit.framework.TestCase;

/**
 */
public class GraphAlignerTest extends TestCase {
  public void testAlign() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "AGCCCAC", "AGCCCC"}, new long[][]{{1, 2}, {1, 3}});
    final Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(1L, 2L), 2, 5, 0, graph));
    final byte[] read = DnaUtils.encodeString("CAAGCCCA");
    checkAlignment(graph, expected, read, 2, new ContigPosition(1, 4, graph));
  }
  GraphAlignment makeAlignment(List<Long> contigs, int start, int end, int score, Graph graph) {
    AlignmentChain current = null;
    for (final long contig: contigs) {
      current = new AlignmentChain(score, new AlignmentSection(contig, start, end), current);
    }
    final AlignmentChain forward = current;
    final AlignmentChain reverse = new AlignmentChain(0, new AlignmentSection(contigs.get(0), start, end), null);
    return new GraphAlignment(forward, reverse, graph);

  }

  private void checkAlignment(Graph graph, Set<GraphAlignment> expected, byte[] read, int readStart, ContigPosition contigStart) {
    checkAlignment(graph, expected, read, readStart, contigStart, 0);
  }
  private void checkAlignment(Graph graph, Set<GraphAlignment> expected, byte[] read, int readStart, ContigPosition contigStart, int mismatches) {
    final GraphAligner aligner = new GraphAligner(graph, new IntegerOrPercentage(mismatches));
    final Set<GraphAlignment> alignments = aligner.align(read, readStart, contigStart);
    assertNotNull("readStart=" + readStart + " contigStart=" + contigStart, alignments);
    assertEquals(expected, alignments);
  }

  public void testGraphAlign2() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "AGCCCGC", "AGCCCC"}, new long[][]{{1, 2}, {1, 3}});
    final byte[] read = DnaUtils.encodeString("CAAGCCCA");
    final GraphAligner aligner = new GraphAligner(graph, new IntegerOrPercentage(0));
    final Collection<GraphAlignment> alignments = aligner.align(read, 2, new ContigPosition(1, 4, graph));
    assertTrue(alignments.isEmpty());
  }

  <T> Set<T> makeSet(T element) {
    final Set<T> s = new HashSet<>();
    s.add(element);
    return s;
  }

  public void testGraphAlign3() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "AGCCCAC", "AGCCCC"}, new long[][]{{1, 2}, {1, 3}});
    //    final List<AlignmentSection> expected = Arrays.asList(new AlignmentSection(1, 2, 5), new AlignmentSection(3, 2, 5));
    final Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(1L, 3L), 2, 5, 0, graph));
    final byte[] read = DnaUtils.encodeString("CAAGCCCC");
    checkAlignment(graph, expected, read, 2, new ContigPosition(1, 4, graph));
  }

  public void testAlignOneReversedContig() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "GTGGGCT", "AGCCCC"}, new long[][]{{1, -2}, {1, 3}});
    //final List<AlignmentSection> expected = Arrays.asList(new AlignmentSection(1, 2, 5), new AlignmentSection(-2, 2, 5));
    final Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(1L, -2L), 2, 5, 0, graph));
    final byte[] read = DnaUtils.encodeString("CAAGCCCA");
    checkAlignment(graph, expected, read, 2, new ContigPosition(1, 4, graph));
  }

  //        |hit here contig -1:position 2
  //      012345
  //      CTTGTT
  // GTGGGCT
  // 0123456
  //  TGGGCTTG  <- read hits at position 6
  //  01234567
  public void testHitToReverse() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "GTGGGCT", "AGCCCC"}, new long[][]{{1, -2}, {1, 3}});
    //final List<AlignmentSection> expected = Arrays.asList(new AlignmentSection(2, 1, 4), new AlignmentSection(-1, 0, 3));
    final Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(2L, -1L), 1, 3, 0, graph));
    final byte[] read = DnaUtils.encodeString("TGGGCTTG");
    checkAlignment(graph, expected, read, 6, new ContigPosition(-1, 2, graph));
  }

  public void testHitAtStartOfRead() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "GTGGGCT", "AGCCCC"}, new long[][]{{1, -2}, {1, 3}});
    //    final List<AlignmentSection> expected = Arrays.asList(new AlignmentSection(2, 1, 6), new AlignmentSection(-1, 2, 3));
    final Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(2L, -1L), 1, 3, 0, graph));
    final byte[] read = DnaUtils.encodeString("TGGGCTTG");
    checkAlignment(graph, expected, read, 0, new ContigPosition(2, 1, graph));
  }
  public void testAllHits() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "GTGGGCT", "AGCCCC"}, new long[][]{{1, -2}, {1, 3}});

    final Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(2L, -1L), 1, 3, 0, graph));
    final byte[] read = DnaUtils.encodeString("TGGGCTTG");
    for (int i = 0; i < read.length; i++) {
      if (i < 6) {
        checkAlignment(graph, expected, read, i, new ContigPosition(2, i + 1, graph));
      }
      if (i > 3) {
        checkAlignment(graph, expected, read, i, new ContigPosition(-1, i - 4, graph));
      }
    }
  }

  public void testAlignFailStart() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "AGCCCAC", "AGCCCC"}, new long[][]{{1, 2}, {1, 3}});
    final byte[] read = DnaUtils.encodeString("GAAGCCCA");
    final GraphAligner aligner = new GraphAligner(graph, new IntegerOrPercentage(0));
    final Collection<GraphAlignment> alignments = aligner.align(read, 2, new ContigPosition(1, 4, graph));
    assertTrue(alignments.isEmpty());
  }

  public void testAlignFailMidway() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"CCCCAC", "CCCCCCCCCCC"}, new long[][]{{1, 2}});
    final byte[] read = DnaUtils.encodeString("CCCCCCCC");
    final GraphAligner aligner = new GraphAligner(graph, new IntegerOrPercentage(0));
    final Collection<GraphAlignment> alignments = aligner.align(read, 2, new ContigPosition(1, 3, graph));
    assertTrue(alignments.isEmpty());
  }
  public void testConnectTo() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "AGCCCAC", "AGCCCC", "GGGG"}, new long[][]{{1, 2}, {1, 3}, {4, 1}});
    //    final List<AlignmentSection> expected = Arrays.asList(new AlignmentSection(1, 2, 5), new AlignmentSection(2, 2, 5));
    final Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(1L, 2L), 2, 5, 0, graph));
    final byte[] read = DnaUtils.encodeString("CAAGCCCA");
    checkAlignment(graph, expected, read, 2, new ContigPosition(1, 4, graph));
  }

  public void testLongRead() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "AGCCCAC", "AGCCCC"}, new long[][]{{1, 2}, {1, 3}});
    final byte[] read = DnaUtils.encodeString("CAAGCCCACCCCC");
    final GraphAligner aligner = new GraphAligner(graph, new IntegerOrPercentage(0));
    final Collection<GraphAlignment> alignments = aligner.align(read, 2, new ContigPosition(1, 3, graph));
    assertTrue(alignments.isEmpty());
  }
  public void testLoop() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"CCCAA", "AAAA", "AATTTT"}, new long[][]{{1, 2}, {2, 3}, {2, 2}});
    //    final List<AlignmentSection> expected = Arrays.asList(new AlignmentSection(1, 1, 2)
    //        , new AlignmentSection(2, 0, 1)
    //        , new AlignmentSection(2, 0, 1)
    //        , new AlignmentSection(2, 0, 3)
    //        , new AlignmentSection(3, 2, 5)
    //    );
    final Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(1L, 2L, 2L, 2L, 3L), 1, 5, 0, graph));
    //                                           AATTTT
    //                                         AAAA
    //                                       AAAA
    //                                     AAAA
    //                                   CCAA
    final byte[] read = DnaUtils.encodeString("CCAAAAAAAATTTT");
    //                                   012345678901234
    checkAlignment(graph, expected, read, 7, new ContigPosition(2, 1, graph));
  }

  public void testMismatchLimit() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AACAAG", "AGCCCAC", "AGCCCC"}, new long[][]{{1, 2}, {1, 3}});
    final Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(1L, 2L), 2, 5, 3, graph));
    //  byte[] read = DnaUtils.encodeString("CAAGCCCA");
    final byte[] read = DnaUtils.encodeString("CTATCTCA");
    checkAlignment(graph, expected, read, 2, new ContigPosition(1, 4, graph), 3);
  }

  public void testPalindrome() {
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"ACGTAA", "AAATTT", "TTGAGA"}, new long[][]{{1, 2}, {1, -2}, {2, 3}, {-2, 3}});
    final Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(1L, 2L, 3L), 0, 5, 0, graph));
    final byte[] read = DnaUtils.encodeString("ACGTAAATTTGAGA");
    checkAlignment(graph, expected, read, 0, new ContigPosition(1, 0, graph), 0);
    final GraphAligner aligner = new GraphAligner(graph, new IntegerOrPercentage(0));
    final Set<GraphAlignment> alignments = aligner.align(read, 0, new ContigPosition(1, 0, graph));

    final byte[] read2 = DnaUtils.encodeString("TCTCAAATTTACGT");
    final Set<GraphAlignment> expected2 = makeSet(makeAlignment(Arrays.asList(-3L, 2L, -1L), 0, 5, 0, graph));
    checkAlignment(graph, expected2, read2, 0, new ContigPosition(-3, 0, graph), 0);
    final Set<GraphAlignment> alignments2 = aligner.align(read2, 0, new ContigPosition(-3, 0, graph));

    assertEquals(1, alignments.size());
    assertEquals(1, alignments2.size());

    /*
    TODO
    GraphAlignment forward = alignments.get(0);
    GraphAlignment reverse = alignments2.get(0);
    int reverseStart = graph.contigLength(reverse.contigs().get(0)) - reverse.startPosition() - 1;
    int reverseEnd = graph.contigLength(reverse.contigs().get(reverse.contigs().size() - 1)) - reverse.endPosition() - 1;
    assertEquals(forward.startPosition(), reverseEnd);
    assertEquals(forward.endPosition(), reverseStart);
    for (int i = 0; i < forward.contigs().size(); i++) {
      long forwardId = forward.contigs().get(i);
      long reverseId = -reverse.contigs().get(reverse.contigs().size() - i - 1);
      assertEquals(forwardId, reverseId);
    }
     */
  }

  public void testSqueezePaths() {
    final Graph graph = GraphMapCliTest.makeGraph(3, new String[]{"ACGTAAT", "AATCTTTG", "TTGAGA"}, new long[][]{{1, 2}, {2, 3}});
    final GraphAligner aligner = new GraphAligner(graph, new IntegerOrPercentage(0));
    final byte[] read = DnaUtils.encodeString("AATCTTTG");
    Set<GraphAlignment> alignments = aligner.align(read, 0, new ContigPosition(1, 4, graph));
    Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(2L), 0, 7, 0, graph));
    assertEquals(expected, alignments);
    alignments = aligner.align(read, 7, new ContigPosition(3, 2, graph));
    assertEquals(expected, alignments);

  }
  public void testHighOverlap() {
    long[][] paths = new long[4][];
    for (int i = 0; i < paths.length; i++) {
      paths[i] = new long[] {i + 1, i + 2};
    }
    final Graph graph = GraphMapCliTest.makeGraph(4, new String[]{"GTAAT", "TAATC", "AATCT", "ATCTT", "TCTTGAGA"}, paths);
    final GraphAligner aligner = new GraphAligner(graph, new IntegerOrPercentage(0));
    final byte[] read = DnaUtils.encodeString("AATCT");
    Set<GraphAlignment> alignments = aligner.align(read, 0, new ContigPosition(1, 2, graph));
    Set<GraphAlignment> expected = makeSet(makeAlignment(Arrays.asList(3L), 0, 4, 0, graph));
    assertEquals(expected, alignments);
    alignments = aligner.align(read, 4, new ContigPosition(5, 2, graph));
    assertEquals(expected, alignments);

  }
}
