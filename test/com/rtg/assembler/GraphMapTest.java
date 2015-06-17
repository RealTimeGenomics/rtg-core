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

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphImplementation;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.NullStreamUtils;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.StandardDeviation;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.store.StoreDirString;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class GraphMapTest extends TestCase {

  private File mDir;

  @Override
  public void setUp() throws IOException {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    Diagnostic.setLogStream();
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  public void testFindPath() {
    final long[][] paths = {{1, 2}, {2, 3}, {1, 2, 3}, {1, 2, 4}, {3, 1, 5}};
    final Graph graph = GraphMapCliTest.makeGraph(3, new String[]{"A", "A", "A", "A", "A"}, paths);
    assertEquals(1, GraphMap.findPath(Arrays.asList(1L, 2L), graph, true));
    assertEquals(-1, GraphMap.findPath(Arrays.asList(-2L, -1L), graph, true));
    assertEquals(-3, GraphMap.findPath(Arrays.asList(-3L, -2L, -1L), graph, true));
    assertEquals(4, GraphMap.findPath(Arrays.asList(1L, 2L, 4L), graph, true));
    assertEquals(0, GraphMap.findPath(Arrays.asList(1L, 2L, 5L), graph, true));
    assertEquals(0, GraphMap.findPath(Arrays.asList(1L, 2L, -3L), graph, true));
    assertEquals(0, GraphMap.findPath(Arrays.asList(4L, 2L, 3L), graph, true));
    assertEquals(0, GraphMap.findPath(Arrays.asList(1L), graph, true));
    assertEquals(0, GraphMap.findPath(Arrays.asList(5L, 1L, 5L), graph, true));
  }

  public void testFindPathPalindrome() {
    final long[][] paths = {{1, 2}, {1, -2}, {2, 3}, {-2, 3}, {1, -2, 3}, {2, 3, 4}};
    final Graph graph = GraphMapCliTest.makeGraph(3, new String[]{"AC", "AT", "AC", "AC"}, paths);
    assertEquals(5, GraphMap.findPath(Arrays.asList(1L, 2L, 3L), graph, true));
    assertEquals(5, GraphMap.findPath(Arrays.asList(1L, -2L, 3L), graph, true));
  }

  public void testFindPathPalindrome2() {
    final long[][] paths = {{1, 2}, {-1, 2}, {2, 3}, {-1, 2, 3}};
    final Graph graph = GraphMapCliTest.makeGraph(3, new String[]{"AT", "AC", "AC"}, paths);
    assertEquals(4, GraphMap.findPath(Arrays.asList(-1L, 2L, 3L), graph, true));
    assertEquals(4, GraphMap.findPath(Arrays.asList(1L, 2L, 3L), graph, true));
  }
  public void testNullRun() throws IOException {
    final GraphImplementation constructGraph = buildGraph(3, "AAAAAAAAAAAAAAAAA");
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    GraphMapTask.run(constructGraph, Collections.<ReadPairSource>emptyList(), params(4, 4, 0), stats);
    TestUtils.containsAll(stats.getStatistics()
        , "0 Mapped without pairing"
        , "0 Too many paths"
        , "0 Too many pairings"
        , "0 No paths"
        , "0 No hits"
        , "0 Mapped in a single contig"
        , "0 Reads mapped across contigs"
        , "0 Equivalent alignments skipped"
    );
  }
  public void testSimplestRun() throws IOException {
    final GraphImplementation constructGraph = buildGraph(3, "AAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGG");
    constructGraph.addPath(new PathArray(1, 2));
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads("AAAAAAAAAAAAGGGGGGGGGGGG")), params(4, 4, 0), stats);
    TestUtils.containsAll(stats.getStatistics()
        , "1 Mapped"
        , "0 Too many paths"
        , "0 No paths"
        , "0 No hits"
        , "0 Mapped in a single contig"
        , "1 Reads mapped across contigs"
        , "25 Equivalent alignments skipped"
    );
    assertEquals("1", constructGraph.pathAttribute(1, GraphKmerAttribute.READ_COUNT));
  }
  public void testSimplestRunRC() throws IOException {
    final GraphImplementation constructGraph = buildGraph(3, "AAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGG");
    constructGraph.addPath(new PathArray(1, 2));

    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads(DnaUtils.reverseComplement("AAAAAAAAAAAAGGGGGGGGGGGG"))), params(4, 4, 0), stats);
    assertEquals("1", constructGraph.pathAttribute(1, GraphKmerAttribute.READ_COUNT));
  }
  public void testNotMappingRun() throws IOException {
    final GraphImplementation constructGraph = buildGraph(3, "AAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGG");
    constructGraph.addPath(new PathArray(1, 2));
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads("AAAAAAAAAAAAGGGGGTATATATATATT"
        , "ATATATTATATATATATATATATAT"
        , "ATATATTATATATATAGGGGGGG"
        , "AAAAATATATATATATATATATAGGGGG"
    )), params(4, 4, 0), stats);
    assertNull(constructGraph.pathAttribute(1, GraphKmerAttribute.READ_COUNT));
    TestUtils.containsAll(stats.getStatistics()
        , "0 Mapped"
        , "0 Too many paths"
        , "4 No paths"
        , "1 No hits"
        , "0 Mapped in a single contig"
        , "0 Reads mapped across contigs"
    );
  }
  public void testTooManyPaths() throws IOException {
    final GraphImplementation constructGraph = buildGraph(3, "AAAAAAAAAAAAAAAAAGGGG", "GGGGGGGGGGGGGGGGGG", "CCCCCCCCCCCCCCCCCC");
    constructGraph.addPath(new PathArray(1, 2));
    constructGraph.addPath(new PathArray(1, -3));
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads("AAAAAAAAAAAAGGGGGGGGGGGGGGGGGG")), params(4, 4, 0), stats);
    TestUtils.containsAll(stats.getStatistics()
        , "0 Mapped"
        , "1 Too many paths"
        , "0 No paths"
        , "0 No hits"
        , "0 Mapped in a single contig"
        , "0 Reads mapped across contigs"
    );
    assertNull(constructGraph.pathAttribute(1, GraphKmerAttribute.READ_COUNT));
    assertNull(constructGraph.pathAttribute(2, GraphKmerAttribute.READ_COUNT));
  }

  public void testMissingLink() throws IOException {
    final GraphImplementation constructGraph = buildGraph(3, "AAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGGGG");
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads("AAAAAAAAAAAAGGGGGGGGGGGGGGGGG")), params(4, 4, 0), stats);
    TestUtils.containsAll(stats.getStatistics()
        , "0 Mapped"
        , "0 Too many paths"
        , "1 No paths"
        , "0 No hits"
        , "0 Mapped in a single contig"
        , "0 Reads mapped across contigs"
    );
  }
  public void testMissedExtension() throws IOException {
    final GraphImplementation constructGraph = buildGraph(3, "CACACACACACAAAA", "AAAAGGGGGGGGGG");
    constructGraph.addPath(new PathArray(1, 2));
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads("CACACACACACAAAAGGCGCGCGCGCGC")), params(4, 4, 0), stats);
    TestUtils.containsAll(stats.getStatistics()
        , "0 Mapped"
        , "0 Too many paths"
        , "1 No paths"
        , "0 No hits"
        , "0 Mapped in a single contig"
        , "0 Reads mapped across contigs"
    );
  }
  GraphMapParams params(int step, int word, int mismatches) {
    return GraphMapParams.builder()
        .stepSize(step)
        .wordSize(word)
        .maxMismatches(new IntegerOrPercentage(mismatches))
        .directory(mDir)
        .create();
  }
  public void testSingleContig() throws IOException {
    final GraphImplementation constructGraph = buildGraph(3, "ATAAAAAAAAAAAAAAA", "CGGGGGGGGGGGGGGG");
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads("TAAAAAAAAAAAAAAA", "TAAAAAAAAAAAAAA", "CCCCCCCCCCCCG")), params(4, 4, 0), stats);
    assertEquals("2", constructGraph.contigAttribute(1, GraphKmerAttribute.READ_COUNT));
    assertEquals("1", constructGraph.contigAttribute(2, GraphKmerAttribute.READ_COUNT));
    TestUtils.containsAll(stats.getStatistics()
        , "3 Mapped"
        , "0 Too many paths"
        , "0 No paths"
        , "0 No hits"
        , "3 Mapped in a single contig"
        , "0 Reads mapped across contigs"
    );
  }
  public void testForwardReverse() throws IOException {
    final GraphImplementation constructGraph = buildGraph(3, "AAACCCTT");
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads("AAACCCTT", "AAGGGTTT")), params(6, 6, 0), stats);
    TestUtils.containsAll(stats.getStatistics()
        , "2 Mapped"
        , "0 Too many paths"
        , "0 No paths"
        , "0 No hits"
        , "2 Mapped in a single contig"
        , "0 Reads mapped across contigs"
    );
    assertEquals("2", constructGraph.contigAttribute(1, GraphKmerAttribute.READ_COUNT));
  }

  public void testNotSimpleSequence() throws IOException {
    final GraphImplementation constructGraph = buildGraph(6, "GCCCTCCGGCCACGCGGAGCCGGCAGTACT"
        , "AGTACTTACTTGCAAGCTGACTGTAGGCACCTTATTTCTAT"
        , "AGTACTTAAAAGTCAACGGTCCACCTGGGACCCAAACGC");
    constructGraph.addPath(new PathArray(1, 2));
    constructGraph.addPath(new PathArray(1, 3));
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    final GraphMapParams params = GraphMapParams.builder()
        .stepSize(6)
        .wordSize(6)
        .maxMismatches(new IntegerOrPercentage(0))
        .directory(mDir)
        .create();
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads(
        "CCACGCGGAGCCGGCAGTACTTACTTGCAAGCTGACTG"
        , "CCACGCGGAGCCGGCAGTACTTAAAAGTCAACGGTCCA")), params, stats);
    TestUtils.containsAll(stats.getStatistics()
        , "2 Mapped"
        , "0 Too many paths"
        , "0 No paths"
        , "0 No hits"
        , "0 Mapped in a single contig"
        , "2 Reads mapped across contigs"
    );
    assertEquals("1", constructGraph.pathAttribute(1, GraphKmerAttribute.READ_COUNT));
    assertEquals("1", constructGraph.pathAttribute(2, GraphKmerAttribute.READ_COUNT));
  }
  public void testNotSimpleSequenceWithErrors() throws IOException {
    final GraphImplementation constructGraph = buildGraph(6, "GCCCTCCGGCCACGCGGAGCCGGCAGTACT"
        , "AGTACTTACTTGCAAGCTGACTGTAGGCACCTTATTTCTAT"
        , "AGTACTTAAAAGTCAACGGTCCACCTGGGACCCAAACGC");
    constructGraph.addPath(new PathArray(1, 2));
    constructGraph.addPath(new PathArray(1, 3));
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    final GraphMapParams params = GraphMapParams.builder()
        .stepSize(6)
        .wordSize(6)
        .maxMismatches(new IntegerOrPercentage(3))
        .directory(mDir)
        .create();
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads(
        "CCACGCGGAGCCGGCAGTACTTTCTTGGAAGCTGACTG"
        , "CCACGCGGAGCCGGCAGTACTTCAATGTCAACGGTCCA")), params, stats);
    TestUtils.containsAll(stats.getStatistics()
        , "2 Mapped"
        , "0 Too many paths"
        , "0 No paths"
        , "0 No hits"
        , "0 Mapped in a single contig"
        , "2 Reads mapped across contigs"
    );
    assertEquals("1", constructGraph.pathAttribute(1, GraphKmerAttribute.READ_COUNT));
    assertEquals("1", constructGraph.pathAttribute(2, GraphKmerAttribute.READ_COUNT));
//    assertEquals(
//        "0/0\tP1\t9\t22" + StringUtils.LS
//        + "1/0\tP2\t9\t22" + StringUtils.LS
//        , alignments.toString());
  }

  ReadPairSource buildReads(String... reads) throws IOException {
    final StringBuilder sb = new StringBuilder();
    final int i = 0;
    for (String read: reads) {
      sb.append(">")
          .append(i)
          .append(StringUtils.LS)
          .append(read)
          .append(StringUtils.LS);
    }
    return new ReadPairSource(ReaderTestUtils.getReaderDnaMemory(sb.toString()));
  }

  private GraphImplementation buildGraph(int overlap, String... contigs) {
    final HashMap<String, String> attr = new HashMap<>();
    attr.put("readCount", "");
    final GraphKmerAttribute constructGraph = new GraphKmerAttribute(overlap, attr, attr);


    for (String contig : contigs) {
      constructGraph.addContig(new ContigString(contig));
    }
    return constructGraph;
  }

  public void testPaired() throws IOException {
    final Map<String, String> readCount = new HashMap<>();
    readCount.put(GraphKmerAttribute.READ_COUNT, "count of reads");
    final MutableGraph graph = GraphMapCliTest.makeGraph(3, new String[]{"ACAACAC", "CGGGGT", "TCCCCTACTACAGCAG"}, new long[][]{{1, 2}, {2, 3}}, readCount, readCount);
    final File tmp = FileHelper.createTempDirectory();
    try {
      final MemoryPrintStream mps =  new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      final SequencesReader left = ReaderTestUtils.getReaderDnaMemory(">left" + LS + "ACAACA" + LS);
      final SequencesReader right = ReaderTestUtils.getReaderDnaMemory(">right" + LS + "CTGCTGTAGTAGGGGA" + LS);
      final ReadPairSource source = new ReadPairSource(left, right);
      source.setMinInsertSize(0);
      source.setMaxInsertSize(20);
      final GraphMap hits = new GraphMap(new GraphIndex(graph, 4, 4), graph, new PairConstraintWriter(mps.printStream()), new PathTracker(new PalindromeTracker(graph)));
      final SimpleThreadPool pool = new SimpleThreadPool(1, "AlignmentIteratorTest", true);
      final AsyncReadSource readSource = new AsyncReadSource(source, "testAlignmentIterator");
      pool.execute(readSource);
      hits.mapReads(readSource, new IntegerOrPercentage(0));
      TestUtils.containsAll(hits.getStatistics().getStatistics(),
          "1 Paired"
          , "0 Too many paths"
          , "0 Too many pairings"
      );
      pool.terminate();
//      assertEquals("0\tP3\t0\t15" + StringUtils.LS, alignments.toString());
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }

  public void testPairedTooManyPathsAndProgress() throws IOException {
    final Map<String, String> readCount = new HashMap<>();
    readCount.put(GraphKmerAttribute.READ_COUNT, "count of reads");
    final MutableGraph graph = GraphMapCliTest.makeGraph(1, new String[]{"CCA", "AAA", "ATT", "ATA"}, new long[][]{{1, 2}, {2, 3}, {2, 2}, {1, 4}, {4, 3}, {2, 4}, {4, 2}}, readCount, readCount);
    final File tmp = FileHelper.createTempDirectory();
    try {
      final MemoryPrintStream mps =  new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      final MemoryPrintStream progress =  new MemoryPrintStream();
      Diagnostic.setProgressStream(progress.printStream());
      final SequencesReader left = ReaderTestUtils.getReaderDnaMemory(">left" + LS + "CCA" + LS);
      final SequencesReader right = ReaderTestUtils.getReaderDnaMemory(">right" + LS + "ATT" + LS);
      final ReadPairSource source = new ReadPairSource(left, right);
      source.setMinInsertSize(5000);
      source.setMaxInsertSize(20000);
      final GraphMap hits = new GraphMap(new GraphIndex(graph, 3, 3), graph, new PairConstraintWriter(NullStreamUtils.getNullPrintStream()), new PathTracker(new PalindromeTracker(graph)));
      final SimpleThreadPool pool = new SimpleThreadPool(1, "AlignmentIteratorTest", true);
      final AsyncReadSource readSource = new AsyncReadSource(source, "testAlignmentIterator");
      pool.execute(readSource);
      hits.mapReads(readSource, new IntegerOrPercentage(0));
      TestUtils.containsAll(hits.getStatistics().getStatistics(),
          "0 Paired"
        , "0 Too many paths"
        , "1 Too many pairings"
        , "0 Equivalent alignments skipped"
          );
      TestUtils.containsAll(progress.toString(),
          "0/1 (0.0%)"
          , "1/1 (100.0%)"
      );
      pool.terminate();
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }

  public void testEquivalentAlignments() throws IOException {
    final GraphImplementation constructGraph = buildGraph(3, "ACGTACGG");
    final StoreDirString graphDir = new StoreDirString();
    GraphWriter.write(constructGraph, graphDir, "foo", Collections.<UUID>emptySet());
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    // Both forward and rc read should have two hits and be aligned once each
    final GraphMapStatistics stats = new GraphMapStatistics(null);
    final GraphMapParams params = params(4, 4, 3);
    GraphMapTask.run(constructGraph, Arrays.asList(buildReads("ACGTACGG", "CCGTACGT")), params, stats);
    assertEquals("2", constructGraph.contigAttribute(1, GraphKmerAttribute.READ_COUNT));
    TestUtils.containsAll(stats.getStatistics()
        , "2 Mapped"
        , "2 Mapped in a single contig"
        , "0 Reads mapped across contigs"
        , "2 Equivalent alignments skipped"
    );
  }

  public void testInsertCalculation() throws IOException {
    final Map<String, String> readCount = new HashMap<>();
    readCount.put(GraphKmerAttribute.READ_COUNT, "count of reads");
    final MutableGraph graph = GraphMapCliTest.makeGraph(3, new String[]{"ACAACAC", "CGGGGT", "TCCCCTACTACAGCAG"}, new long[][]{{1, 2}, {2, 3}}, readCount, readCount);
    final File tmp = FileHelper.createTempDirectory();
    try {
      final MemoryPrintStream mps =  new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      final SequencesReader left = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("TCCCCTACT", "CCTACTA"));
      final SequencesReader right = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("CTGCTGTAGTA", "CTGCTGTA"));
      final ReadPairSource source = new ReadPairSource(left, right);
      final GraphMap hits = new GraphMap(new GraphIndex(graph, 4, 4), graph, new PairConstraintWriter(mps.printStream()), new PathTracker(new PalindromeTracker(graph)));
      final SimpleThreadPool pool = new SimpleThreadPool(1, "AlignmentIteratorTest", true);
      final AsyncReadSource readSource = new AsyncReadSource(source, "testAlignmentIterator");
      pool.execute(readSource);
      pool.terminate();
      final StandardDeviation dev = hits.calculateInserts(readSource, new IntegerOrPercentage(0));
      assertEquals(-2.0, dev.mean());
      assertEquals(1.0, dev.standardDeviation());
    } finally {
      FileHelper.deleteAll(tmp);
    }
  }
}
