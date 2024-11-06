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
    final SimpleThreadPool pool = new SimpleThreadPool(1, "AlignmentIteratorTest", true);
    final AsyncReadSource readSource = new AsyncReadSource(readPairSource, "testAlignmentIterator");
    pool.execute(readSource);
    final AlignmentIterator iterator = new AlignmentIterator(readSource, graph, graphAligner, graphIndex, graphMapStatistics);
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
