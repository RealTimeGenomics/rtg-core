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

package com.rtg.assembler.graph.implementation;

import java.util.HashMap;
import java.util.Map;

import com.rtg.assembler.graph.Contig;

import junit.framework.TestCase;

/**
 */
public class GraphKmerAttributeTest extends TestCase {

  static GraphKmerAttribute graph() {
    final String[] contigNt = {"ACGT", "CGTA", "GTAC", "TACG", "AA", "CC", "GG", "TT"};
    final Contig[] contigs = new Contig[contigNt.length];
    for (int i = 0; i < contigNt.length; i++) {
      contigs[i] = new ContigString(contigNt[i]);
    }

    final long[] p1 = {+1, -2};
    final long[] p2 = {+2, -3, +4};
    final long[] p3 = {-3, +5, +6};
    final long[] p4 = {+3, +7};
    final long[][] paths = {p1, p2, p3, p4};

    final Map<String, String> contigAttributes = new HashMap<>();
    contigAttributes.put("foo", "foo comment");
    contigAttributes.put("bar", "bar comment");

    final Map<String, String> pathAttributes = new HashMap<>();
    pathAttributes.put("boo", "comment boo");
    pathAttributes.put("oob", "comment oob");

    final GraphKmerAttribute graph = new GraphKmerAttribute(0, contigAttributes, pathAttributes);
    GraphImplementationTest.buildGraph(graph, contigs, paths);
    return graph;
  }

  public void test() {
    final GraphKmerAttribute graph = graph();
    assertEquals(8, graph.numberContigs());

    graph.setKmerFreq(1, 42);
    graph.setContigAttribute(2, "kMerFreq", "123");
    graph.setContigAttribute(2, "foo", "baz");
    graph.setKmerFreq(-8, 101);

    assertEquals(42, graph.kmerFreq(1));
    assertEquals("42", graph.contigAttribute(-1, "kMerFreq"));

    assertEquals(123, graph.kmerFreq(-2));
    assertEquals("123", graph.contigAttribute(2, "kMerFreq"));
    assertEquals("baz", graph.contigAttribute(2, "foo"));
    assertEquals("baz", graph.contigAttribute(-2, "foo"));
    assertEquals(null, graph.contigAttribute(2, "bar"));

    assertEquals(0, graph.kmerFreq(-3));
    assertEquals("0", graph.contigAttribute(3, "kMerFreq"));

    assertEquals(101, graph.kmerFreq(8));
    assertEquals("101", graph.contigAttribute(-8, "kMerFreq"));


  }
}
