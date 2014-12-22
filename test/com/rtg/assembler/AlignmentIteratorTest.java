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
import java.util.Arrays;

import com.rtg.assembler.graph.Graph;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class AlignmentIteratorTest extends TestCase {
  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void testAlignmentIterator() throws IOException {
    final ReadPairSource readPairSource = new ReadPairSource(
        ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("AGTACGTAGGGATCC", "ACGCAGCGAT", "AAAA"))
        , ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("AACCTTATACA", "AAACCTTAGG", "AGGG")));
    final Graph graph = GraphMapCliTest.makeGraph(2, new String[]{"AAACACCAAGTACGTAGGGATCCAC", "ACACCAACCCTTATACAGTG"}, new long[][]{{1, 2}});
    final GraphMapStatistics graphMapStatistics = new GraphMapStatistics(null);
    final IntegerOrPercentage mismatches = new IntegerOrPercentage(3);
    final GraphAligner graphAligner = new GraphAligner(graph, mismatches, new GraphTraversions(graph));
    final GraphIndex graphIndex = new GraphIndex(graph, 5, 5);
    SimpleThreadPool pool = new SimpleThreadPool(1, "AlignmentIteratorTest", true);
    final AsyncReadSource readSource = new AsyncReadSource(readPairSource, "testAlignmentIterator");
    pool.execute(readSource);
    AlignmentIterator iterator = new AlignmentIterator(readSource, graph, graphAligner, graphIndex, graphMapStatistics);
    assertTrue(iterator.hasNext());
    AlignmentIterator.ReadAlignment next = iterator.next();
    assertEquals(0, next.mId);
    assertEquals(2, next.mFragments.size());
    assertTrue(next.mFragments.get(0).toString(), next.mFragments.get(0).contains(new GraphAlignment(8, 22, Arrays.asList(1L), 0, graph)));
    assertTrue(next.mFragments.get(1).toString(), next.mFragments.get(1).contains(new GraphAlignment(6, 16, Arrays.asList(2L), 1, graph)));

    assertTrue(iterator.hasNext());
    next = iterator.next();

    assertEquals(1, next.mId);
    assertEquals(2, next.mFragments.size());
    assertTrue(next.mFragments.get(0).toString(), next.mFragments.get(0).isEmpty());
    assertTrue(next.mFragments.get(1).toString(), next.mFragments.get(1).isEmpty());
    assertTrue(iterator.hasNext());
    next = iterator.next();

    assertEquals(2, next.mId);
    assertEquals(2, next.mFragments.size());
    assertTrue(next.mFragments.get(0).toString(), next.mFragments.get(0).isEmpty());
    assertTrue(next.mFragments.get(1).toString(), next.mFragments.get(1).isEmpty());

    assertFalse(iterator.hasNext());
    pool.terminate();

    TestUtils.containsAll(graphMapStatistics.getStatistics()
        , "4 No hits"
        , "1 Equivalent alignments skipped"
    );
  }
}
