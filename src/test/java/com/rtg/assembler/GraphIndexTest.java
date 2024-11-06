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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphImplementation;
import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 */
public class GraphIndexTest extends TestCase {
  static final String[] CONTIGS = {
     //12345678|2345678|2345678|2345678|23456 wordsize 8
      "ACCCGGGTTGAGAGGTGTGGGTGGTGACGTACGTACGT"
       //12345678|2345678|2345678|2345678|23456 wordsize 8
      , "TCTCGCGTTAAGCGGAGTGAGTGATGACGTACGTACGT"

  };

  public void testBuild() throws IOException {
    final MutableGraph graph = makeGraph(CONTIGS);
    final TreeMap<Long, Long> positionDecoder = ContigPosition.buildDecoder(graph);
    final GraphIndex index = buildIndex(graph);
    List<Long> hits = index.hitContigs(Long.parseLong("01112223", 4));
    assertEquals(1, hits.size());
    assertEquals(1, new ContigPosition(hits.get(0), graph, positionDecoder).mContigId);
    assertEquals(7, new ContigPosition(hits.get(0), graph, positionDecoder).mPosition);

    hits = index.hitContigs(Long.parseLong("32012301", 4));
    assertEquals(2, hits.size());

    assertEquals(1, new ContigPosition(hits.get(0), graph, positionDecoder).mContigId);
    assertEquals(31, new ContigPosition(hits.get(0), graph, positionDecoder).mPosition);

    assertEquals(2, new ContigPosition(hits.get(1), graph, positionDecoder).mContigId);
    assertEquals(31, new ContigPosition(hits.get(1), graph, positionDecoder).mPosition);
  }
  public void testNumberOfBits() throws IOException {
    MutableGraph graph = makeGraph(new String[] {"AAAACCCCGGGGTTTT"});
    TreeMap<Long, Long> positionDecoder = ContigPosition.buildDecoder(graph);
    GraphIndex index = buildIndex(graph);
    List<Long> longs = index.hitContigs(Long.parseLong("22223333", 4));
    assertEquals(1, longs.size());
    ContigPosition contigPosition = new ContigPosition(longs.get(0), graph, positionDecoder);
    assertEquals(1, contigPosition.mContigId);
    assertEquals(15, contigPosition.mPosition);

    graph = makeGraph(new String[] {"AAAACCCCGGGGTTTTG"});
    positionDecoder = ContigPosition.buildDecoder(graph);
    index = buildIndex(graph);
    longs = index.hitContigs(Long.parseLong("22223333", 4));
    assertEquals(1, longs.size());
    contigPosition = new ContigPosition(longs.get(0), graph, positionDecoder);
    assertEquals(1, contigPosition.mContigId);
    assertEquals(15, contigPosition.mPosition);
  }

  public static GraphIndex buildIndex(Graph graph) {
    return new GraphIndex(graph, 8, 8);
  }

  static MutableGraph makeGraph(String[] contigs) {
    final Map<String, String> empty = new HashMap<>();
    final MutableGraph graph = new GraphImplementation(2, empty, empty);
    for (String contig : contigs) {
      graph.addContig(new ContigString(contig));
    }
    return graph;
  }

  public void testSearch() throws IOException {
    final MutableGraph graph = GraphIndexTest.makeGraph(GraphIndexTest.CONTIGS);
    final GraphIndex index = GraphIndexTest.buildIndex(graph);
    //            |       |       |
    String read = "TGGGTGGTGACGTACGTACGT";
    List<List<ContigPosition>> contigs = index.hits(DnaUtils.encodeString(read), graph, index.getSearchFunction());
    final Map<Integer, List<ContigPosition>> expected = new HashMap<>();
    expected.put(14, Arrays.asList(new ContigPosition(1L, 31, graph), new ContigPosition(2L, 31, graph)));
    check(expected, contigs, read.length());

    expected.clear();
//    expected.put(13, Arrays.asList(-1L, -2L));
    expected.put(13, Arrays.asList(new ContigPosition(-1L, 13, graph), new ContigPosition(-2L, 13, graph)));
    contigs = index.hits(DnaUtils.encodeString(DnaUtils.reverseComplement(read)), graph, index.getSearchFunction());
    check(expected, contigs, read.length());

    //     |       |    *  |
    read = "TGGGTGGTGACGAACGTACGT";
    contigs = index.hits(DnaUtils.encodeString(read), graph, index.getSearchFunction());
    expected.clear();
    check(expected, contigs, read.length());

    //      |       |       |
    read = "GTGGGTGGTGACGTACGTACGT";
    contigs = index.hits(DnaUtils.encodeString(read), graph, index.getSearchFunction());
    expected.clear();
//    expected.put(7, Arrays.asList(1L));
    expected.put(7, Arrays.asList(new ContigPosition(1L, 23, graph)));
//    expected.put(15, Arrays.asList(1L, 2L));
    expected.put(15, Arrays.asList(new ContigPosition(1L, 31, graph), new ContigPosition(2L, 31, graph)));
    check(expected, contigs, read.length());
  }
  void check(Map<Integer, List<ContigPosition>> expected, List<List<ContigPosition>> actual, int length) {
    assertEquals(length, actual.size());
    for (int i = 0; i < actual.size(); ++i) {
      List<ContigPosition> currentExpected = expected.get(i);
      if (currentExpected == null) {
        currentExpected = Collections.emptyList();
      }
      assertEquals("At: " + i + " actually: " + actual, currentExpected, actual.get(i));
    }

  }


}
