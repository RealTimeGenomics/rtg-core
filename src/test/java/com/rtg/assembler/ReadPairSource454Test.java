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

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class ReadPairSource454Test extends TestCase {
  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }
  public void testAligner() throws IOException {
    final SequencesReader reader1 = ReaderTestUtils.getReaderDnaMemory(lines(">1", "AAACCCTTG"));
    final GraphKmerAttribute graph = GraphMapCliTest.makeGraph(2, new String[]{}, new long[][]{});
    assertTrue(new ReadPairSource454(reader1).aligner(graph, new IntegerOrPercentage(2), new GraphTraversions(graph)) instanceof GraphFlowAligner);
  }
  static String lines(String... lines) {
    final StringBuilder sb = new StringBuilder();
    for (String line : lines) {
      sb.append(line);
      sb.append(LS);
    }
    return sb.toString();
  }
  public void testReads() throws IOException {
    final SequencesReader reader1 = ReaderTestUtils.getReaderDnaMemory(lines(">1", "AAACCCTTG"));
    final SequencesReader reader2 = ReaderTestUtils.getReaderDnaMemory(lines(">1", "ATACGCTTG"));
    final ReadPairSource source = new ReadPairSource454(reader1, reader2);
    final List<byte[]> fragments = source.nextFragments();
    assertNotNull(fragments);
    ReadPairSourceTest.listEquals(Arrays.asList(DnaUtils.encodeString("AAACCCTTG")
        , DnaUtils.encodeString("CAAGCGTAT")
    ), fragments);
    assertNull(source.nextFragments());
  }
}
