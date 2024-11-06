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

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class GraphMapTaskTest extends TestCase {
  public void testReadHitsTask() throws IOException {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File reads = new File(tmpDir, "reads");
      final File graph = new File(tmpDir, "graph");
      assertTrue(graph.mkdir());
      final File output = new File(tmpDir, "output");
      assertTrue(output.mkdir());
      final MutableGraph g = GraphMapCliTest.makeGraph(4, new String[]{"ACGTTTT", "TTTTCC", "TTTTAA"}, new long[][]{{1, 2}, {1, 3}});
      ReaderTestUtils.createPairedReaderDNA(ReadPairSourceTest.LEFT_SEQUENCE, ReadPairSourceTest.RIGHT_SEQUENCE, reads, new SdfId());
      final GraphMapParams params = GraphMapParams.builder()
          .reads(Arrays.asList(reads))
          .graph(g)
          .directory(output)
          .wordSize(4)
          .stepSize(4)
          .create();
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      try {
      new GraphMapTask(params, out.printStream()).run();
      assertTrue(new File(output, "contig.1.fa").exists());
      assertTrue(new File(output, "path.1.tsv").exists());
      assertTrue(new File(output, "header.tsv").exists());

      TestUtils.containsAll(mps.toString()
          , "1 Too many paths"
          , "3 No paths"
      );
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }
  public void testReadHitsTask454() throws IOException {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final File reads = new File(tmpDir, "reads");
      final File graph = new File(tmpDir, "graph");
      assertTrue(graph.mkdir());
      final File output = new File(tmpDir, "output");
      assertTrue(output.mkdir());
      final MutableGraph g = GraphMapCliTest.makeGraph(4, new String[]{"ACGTTTT", "TTTTCC", "TTTTAA"}, new long[][]{{1, 2}, {1, 3}});
      ReaderTestUtils.createPairedReaderDNA(ReadPairSourceTest.LEFT_SEQUENCE, ReadPairSourceTest.RIGHT_SEQUENCE, reads, new SdfId());
      final GraphMapParams params = GraphMapParams.builder()
          .reads454(Arrays.asList(reads))
          .graph(g)
          .directory(output)
          .wordSize(4)
          .stepSize(4)
          .create();
      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      try {
        new GraphMapTask(params, TestUtils.getNullOutputStream()).run();
        assertTrue(new File(output, "contig.1.fa").exists());
        assertTrue(new File(output, "path.1.tsv").exists());
        assertTrue(new File(output, "header.tsv").exists());

        TestUtils.containsAll(mps.toString()
            , "1 Too many paths"
            , "3 No paths"
        );
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }

  public void testReadHitsTaskCalculateInsert() throws IOException {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final Map<String, String> readCount = new HashMap<>();
      readCount.put(GraphKmerAttribute.READ_COUNT, "count of reads");
      final MutableGraph g = GraphMapCliTest.makeGraph(4, new String[]{"ACAACAC", "CGGGGT", "TCCCCTACTACAGCAG"}, new long[][]{{1, 2}, {2, 3}}, readCount, readCount);
      final MemoryPrintStream mps =  new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      final File reads = new File(tmpDir, "reads");
      final File reads2 = new File(tmpDir, "reads2");
      final File reads3 = new File(tmpDir, "reads3");
      final File graph = new File(tmpDir, "graph");
      assertTrue(graph.mkdir());
      final File output = new File(tmpDir, "output");
      assertTrue(output.mkdir());
      ReaderTestUtils.createPairedReaderDNA(ReaderTestUtils.fasta("TCCCCTACT", "CCTACTA"), ReaderTestUtils.fasta("CTGCTGTAGTA", "CTGCTGTA"), reads, new SdfId());
      ReaderTestUtils.createPairedReaderDNA(ReaderTestUtils.fasta("TCCCCTACT", "CCTACTA"), ReaderTestUtils.fasta("CTGCTGTAGTA", "CTGCTGTA"), reads2, new SdfId());
      // make single end long enough that it will time out if blocked
      final String[] singleEndStrings = new String[1049];
      for (int i = 0; i < 1049; ++i) {
        singleEndStrings[i] = "TCCCCTACT";
      }
      ReaderTestUtils.getDNADir(ReaderTestUtils.fasta(singleEndStrings), reads3);
      final GraphMapParams params = GraphMapParams.builder()
          .reads(Arrays.asList(reads, reads2, reads3))
          .graph(g)
          .directory(output)
          .wordSize(4)
          .stepSize(4)
          .create();
      try {
        new GraphMapTask(params, TestUtils.getNullOutputStream()).run();
        assertTrue(new File(output, "contig.1.fa").exists());
        assertTrue(new File(output, "path.1.tsv").exists());
        assertTrue(new File(output, "header.tsv").exists());

        TestUtils.containsAll(mps.toString()
            , "4 Paired end reads"
            , "4 Successfully paired"
            , "Min Insert: -8"
            , "Max Insert: 4"
        );
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }
}
