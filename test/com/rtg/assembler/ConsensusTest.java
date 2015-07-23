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
import com.rtg.assembler.graph.io.GraphReader;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.NullStreamUtils;
import com.rtg.util.store.StoreDirString;

/**
 */
public class ConsensusTest extends AbstractCliTest {

  public void testHelp() {
    checkHelp("rtg consensus"
        , "input graph directory"
        , "-o,", "--output", "output directory"
        , "--consensus-reads", "number of reads necessary to form consensus along a path"
        , "-k,", "--kmer-size", "size of kmer used to build the graph"
    );
  }

  @Override
  protected AbstractCli getCli() {
    return new Consensus();
  }

  public void testBuildConsensus() throws IOException {
    final StoreDirString input = new StoreDirString();
    final StoreDirString output = new StoreDirString();
    final Map<String, String> attr = new HashMap<>();
    attr.put(GraphKmerAttribute.READ_COUNT, "foo");

    final MutableGraph graph = GraphMapCliTest.makeGraph(1, new String[]{"ACGT", "TCGA", "ACCCC", "CGCGCC", "GGGT", "CTTTT"}
        , new long[][]{{1, 2, 3, 4}, {5, 2, 3, 6}}
        , attr
        , attr
    );
    graph.setPathAttribute(1, "readCount", "3");
    graph.setPathAttribute(2, "readCount", "3");
    final Set<UUID> set = new HashSet<>();
    GraphWriter.write(graph, input, "written", set);

    Consensus.writeConsensus(2, 3, input, output);
    final MutableGraph result = (MutableGraph) GraphReader.read(output);
//    System.err.println(graphString(result));
    for (long i = 1; i <= 6; i++) {
      assertTrue(result.contigDeleted(i));
    }
    final Graph compact = GraphSorter.sortedGraph(result.compact());
    assertEquals(2L, compact.numberContigs());
    assertEquals(graphString(result), "AAAAGGGGTCGACCC", ContigString.contigSequenceString(compact.contig(1)));
    assertEquals(graphString(result), "ACGTCGACCCCGCGCC", ContigString.contigSequenceString(compact.contig(2)));

  }

  public void testNxGraph() {
    final MutableGraph graph = GraphMapCliTest.makeGraph(1, new String[]{"ACGTTCGAACCC", "CGCC", "GGGT", "CTTT"}
        , new long[][]{}
    );
    final int[] h = Consensus.nxGraph(graph);
    assertEquals(4, h[100]);
    assertEquals(4, h[50]);
    assertEquals(12, h[49]);
  }
  public void testNxGraph2() {
    final MutableGraph graph = GraphMapCliTest.makeGraph(1, new String[]{"ACGTTCGAACCC", "CGCCGGGT", "CTTT"}
        , new long[][]{}
    );
    final int[] h = Consensus.nxGraph(graph);
    assertEquals(4, h[100]);
    assertEquals(4, h[84]);
    assertEquals(8, h[83]);
    assertEquals(8, h[51]);
    assertEquals(8, h[50]);
    assertEquals(12, h[49]);
  }

  boolean reversed(long permutation, long contig) {
    return (permutation & 1 << (contig - 1)) != 0;
  }

  long factorial(long n) {
    long i = n;
    long fact = 1;
    while (i > 1) {
      fact = fact * i;
      i--;
    }
    return fact;
  }
  public void testFactorial() {
    assertEquals(1, factorial(1));
    assertEquals(2, factorial(2));
    assertEquals(6, factorial(3));
    assertEquals(24, factorial(4));
  }

  private GraphKmerAttribute permutation(GraphKmerAttribute original, long permutation) {
    final GraphSorter.NegativeMap translate = new GraphSorter.NegativeMap();
    final List<Long> unused = new ArrayList<>();
    final List<Long> order = new ArrayList<>();
    for (long j = 1; j <= original.numberContigs(); j++) {
      unused.add(j);
    }

    long remainder = permutation;
    while (unused.size() > 1) {
      final long factorial = factorial(unused.size() - 1);
      final long currentContig = remainder / factorial;
      final Long key = unused.get((int) currentContig);
      translate.put(key, (long) unused.size());
      order.add(0, key);
      unused.remove((int) currentContig);
      remainder = remainder % factorial;
    }
    order.add(0, unused.get(0));
    translate.put(unused.get(0), 1L);

    final GraphKmerAttribute current = new GraphKmerAttribute(original.contigOverlap(), original.contigAttributes(), original.pathAttributes());
    for (int i = 0; i < original.numberContigs(); i++) {
      current.addContig(original.contig(order.get(i)));
    }
    for (long path = 1; path <= original.numberPaths(); path++) {
      final long[] contigs = new long[original.pathLength(path)];
      for (int contig = 0; contig < original.pathLength(path); contig++) {
        contigs[contig] = translate.get(original.pathContig(path, contig));
      }
      current.addPath(new PathArray(contigs));
    }
    return current;
  }
  private GraphKmerAttribute complement(GraphKmerAttribute original, long i) {
    final GraphKmerAttribute current = new GraphKmerAttribute(original.contigOverlap(), original.contigAttributes(), original.pathAttributes());
    for (long contig = 1; contig <= original.numberContigs(); contig++) {
      if (reversed(i, contig)) {
        current.addContig(original.contig(-contig));
      } else {
        current.addContig(original.contig(contig));
      }
    }
    for (long path = 1; path <= original.numberPaths(); path++) {
      final long[] contigs = new long[original.pathLength(path)];
      for (int contig = 0; contig < original.pathLength(path); contig++) {
        final long contigId = original.pathContig(path, contig);
        if (reversed(i, contigId)) {
          contigs[contig] = -contigId;
        } else {
          contigs[contig] = contigId;
        }
      }
      current.addPath(new PathArray(contigs));
    }
    return current;
  }

  public void testPalindromesEndToEnd() throws IOException {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{"TAACGAACCGG", "ACCGGT", "CCGGTCCAGTA"}, new long[][]{{1, 2}, {1, -2}, {2, 3}, {-2, 3}});
    graph.addContigAttribute(GraphKmerAttribute.READ_COUNT, "r");
    graph.addPathAttribute(GraphKmerAttribute.READ_COUNT, "r");
    final double pow = Math.pow(2, graph.numberContigs() + 1);
    final long fact = factorial(graph.numberContigs());
    for (long i = 0; i < pow; i++) {
      for (long j = 0; j < fact; j++) {
        final GraphKmerAttribute current = permutation(complement(graph, i), j);
        final SequencesReader reads = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("TAACGAACCGGTCCAGTA", "TACTGGACCGGTTCGTTA"));
        final PathTracker pathTracker = new PathTracker(new PalindromeTracker(current));
        final GraphMap mapper = new GraphMap(new GraphIndex(current, 5, 5), current, new PairConstraintWriter(NullStreamUtils.getNullPrintStream()), pathTracker);
        final AsyncReadPool readPool = new AsyncReadPool("BAR", Collections.singletonList(new ReadPairSource(reads)));
        mapper.mapReads(readPool.sources().get(0), new IntegerOrPercentage(0));
//      mapper.getStatistics().printStatistics(System.err);
        PathTracker.apply(PathTracker.merge(Collections.singletonList(pathTracker)), current);


        Consensus.buildConsensus(current.contigOverlap(), 1, current);
//      printGraph(current);
        final MutableGraph compact = current.compact();
//      printGraph(compact);
        assertEquals(1, compact.numberContigs());
        assertEquals(0, compact.numberPaths());
        readPool.terminate();
      }
    }

  }
  public void testPalindromes() {
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(5, new String[]{"AATATAATACTG", "GGGTAACG", "ACCGGT", "CCGGTTCGTTA", "CCGGTCCAGTA"}, new long[][] {{5, -1}, {4, -2}, {3, 4}, {-3, 4}, {3, 5}, {-3, 5}, {-5, 3, 4}});
    graph.setPathAttribute(1, GraphKmerAttribute.READ_COUNT, "25");
    graph.setPathAttribute(2, GraphKmerAttribute.READ_COUNT, "25");
    graph.setPathAttribute(4, GraphKmerAttribute.READ_COUNT, "5");
    graph.setPathAttribute(6, GraphKmerAttribute.READ_COUNT, "2");
    graph.setPathAttribute(7, GraphKmerAttribute.READ_COUNT, "34");
//      printGraph(current);
        FilterPaths.improveSingle(graph, 3);
        FilterPaths.improveMultiple(graph, 1);
        Consensus.buildConsensus(graph.contigOverlap(), 5, graph);
//      printGraph(current);
//    System.err.print(graphString(graph));
        final MutableGraph compact = graph.compact();
//      printGraph(compact);
//    System.err.print(graphString(compact));
        assertEquals(1, compact.numberContigs());
        assertEquals(0, compact.numberPaths());
  }
  private String graphString(Graph compact) throws IOException {
    final StoreDirString dir2 = new StoreDirString();
    GraphWriter.writeWithDeleted(compact, dir2, "foo", Collections.<UUID>emptySet());
    return dir2.toString();
  }
}
