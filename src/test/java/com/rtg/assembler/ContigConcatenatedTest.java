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

import java.util.ArrayList;
import java.util.List;

import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;

import junit.framework.TestCase;

/**
 */
public class ContigConcatenatedTest extends TestCase {
  public void test() {
    final GraphKmerAttribute graph = new GraphKmerAttribute(2);
    final List<Long> concatenate = new ArrayList<>();
    concatenate.add(graph.addContig(new ContigString("GACGT")));
    concatenate.add(graph.addContig(new ContigString("GTCCGT")));
    concatenate.add(graph.addContig(new ContigString("GTAC")));
    assertEquals("GACGTCCGTAC", ContigString.contigSequenceString(new ContigConcatenated(concatenate, graph, 3)));
  }

  public void testDifferentK() {
    final GraphKmerAttribute graph = new GraphKmerAttribute(2);
    final List<Long> concatenate = new ArrayList<>();
    concatenate.add(graph.addContig(new ContigString("GACGTC")));
    concatenate.add(graph.addContig(new ContigString("GTCCGTA")));
    concatenate.add(graph.addContig(new ContigString("GTAC")));
    assertEquals("GACGTCCGTAC", ContigString.contigSequenceString(new ContigConcatenated(concatenate, graph, 4)));
  }
  public void testSingle() {
    final GraphKmerAttribute graph = new GraphKmerAttribute(2);
    final List<Long> concatenate = new ArrayList<>();
    concatenate.add(graph.addContig(new ContigString("GACGTC")));
    assertEquals("GACGTC", ContigString.contigSequenceString(new ContigConcatenated(concatenate, graph, 4)));
  }
}
