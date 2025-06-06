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
