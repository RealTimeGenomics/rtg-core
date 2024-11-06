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

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

/**
 */
public class GraphToPlotTest extends AbstractCliTest {

  public void testHelp() {
    checkHelp("rtg graph2plot"
        , "Produces graphs of the contigs contained in the specified graph directory"
        , "input graph directory"
        , "-o,", "--output"
        , "-w", "--width", "maximum distance from the initial node within the .dot"
        , "-s", "--start", "produce a .dot file for nodes around this one"
    );
  }

  public void test() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      TestUtils.containsAll(checkHandleFlagsErr(),
        "Usage: rtg graph2plot [OPTION]... -o DIR -s INT DIR"
        , "You must provide values for -o DIR -s INT DIR"
      );
      assertTrue(tmpDir.list().length == 0);

      final File inFile = new File(tmpDir, "in");
      final File outFile = new File(tmpDir, "out");
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr(inFile.toString(), "-o", outFile.toString(), "-s", "4"), "Input file should be a directory in the RTG graph file format");
      assertTrue(tmpDir.list().length == 0);

      assertTrue(inFile.mkdir());
      assertTrue(outFile.createNewFile());
      TestUtils.containsAllUnwrapped(checkHandleFlagsErr(inFile.toString(), "-o", outFile.toString(), "-s", "4"), "The directory \"" + outFile + "\" already exists.");

      assertTrue(outFile.delete());
      checkHandleFlagsOut(inFile.toString(), "-o", outFile.toString(), "-w", "5", "-s", "2");
      FileHelper.deleteAll(outFile);
      //new GraphToPlot().mainInit(new String[] {inFile.toString(), outFile.toString()}, out.printStream(), err.printStream());
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
  static final String EXPECTED_DOTFILE = ""
      + "digraph contigGraph {" + LS
      + "graph [rankdir=LR, ratio=fill]" + LS
      + "node0 [label=\"+1\\n4   2\"];" + LS
      + "node1 [label=\"-2\\n5   3\"];" + LS
      + "node2 [label=\"+3\\n6   4\", color=red, shape=box];" + LS
      + "node3 [label=\"+4\\n7   5\"];" + LS
      + "node4 [label=\"+5\\n8   6\"];" + LS
      + "node5 [label=\"+6\\n9   7\"];" + LS
      + "node6 [label=\"+7\\n10   8\", style=dashed];" + LS
      + "node7 [label=\"+9\\n12   10\", style=dashed];" + LS
      + "node1 -> node4[];" + LS
      + "node2 -> node4[];" + LS
      + "node1 -> node2[];" + LS
      + "node5 -> node7[];" + LS
      + "node2 -> node5[];" + LS
      + "node4 -> node7[];" + LS
      + "node6 -> node5[];" + LS
      + "node0 -> node1[label=\"2\"];" + LS
      + "node3 -> node2[];" + LS
      + "}" + LS
      + "";
  static final String EXPECTED_LONG_PATHS = ""
      + "{simpleCount: 0, overlapCount: 0, contigs: [-2, 3, 5]}," + LS
      + "{simpleCount: 0, overlapCount: 0, contigs: [4, 3, 6]}" + LS
      + "";

  public void testDotFile() throws IOException {
    final MemoryPrintStream graphOut = new MemoryPrintStream();
    final MemoryPrintStream linksOut = new MemoryPrintStream();
    final Map<String, String> contigAttr = new HashMap<>();
    final Map<String, String> pathAttr = new HashMap<>();
    contigAttr.put(GraphKmerAttribute.READ_COUNT, "desc");
    pathAttr.put(GraphKmerAttribute.READ_COUNT, "desc");
    final GraphKmerAttribute graph = new GraphKmerAttribute(2, contigAttr, pathAttr);
    for (int i = 0; i < 10; ++i) {
      final int length = i + 4;
      contig(graph, length, length * (i + 1));
    }
    graph.addPath(new PathArray(1, -2));
    graph.addPath(new PathArray(-2, 3, 5));
    graph.addPath(new PathArray(-2, 5));
    graph.addPath(new PathArray(4, 3, 6));
    graph.addPath(new PathArray(7, 6));
    graph.addPath(new PathArray(7, 8));
    graph.addPath(new PathArray(5, 9));
    graph.addPath(new PathArray(6, 9));
    graph.addPath(new PathArray(9, 10));
    graph.setPathAttribute(1, GraphKmerAttribute.READ_COUNT, "2");
    graph.setContigAttribute(1, GraphKmerAttribute.READ_COUNT, "3");
    GraphToPlot.writeDotFile(graph, graphOut.printStream(), linksOut.printStream(), 2, 3L);
    mNano.check("graph-to-plot-dot.txt", graphOut.toString());
    //assertEquals(EXPECTED_DOTFILE, graphOut.toString());
    assertEquals(EXPECTED_LONG_PATHS, linksOut.toString());

  }
  void contig(GraphKmerAttribute graph, int length, int kmer) {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < length; ++i) {
      sb.append('A');
    }
    final long id = graph.addContig(new ContigString(sb.toString()));
    graph.setKmerFreq(id, kmer);


  }
  public void testOverlappingCount() {
    final MutableGraph graph = makeGraph(
        new String[]{"AAAAA", "CCCC", "GGGG", "TTTT", "AACC", "AAGG", "AATT"}
        , new long[][]{
             {1, 2}     // 1
        ,       {2, 3}  // 2
        , {4, 1, 2}     // 3
        , {4, 1, 2, 3}  // 4
        , {4, 1}        // 5
        , {6, 1}        // 6
        , {4, 1, 5}     // 7
        , {6, 1, 2, 3}  // 8
    });
    graph.setPathAttribute(1, GraphKmerAttribute.READ_COUNT, "1");
    graph.setPathAttribute(2, GraphKmerAttribute.READ_COUNT, "2");
    graph.setPathAttribute(3, GraphKmerAttribute.READ_COUNT, "3");
    graph.setPathAttribute(4, GraphKmerAttribute.READ_COUNT, "4");
    graph.setPathAttribute(5, GraphKmerAttribute.READ_COUNT, "5");
    graph.setPathAttribute(6, GraphKmerAttribute.READ_COUNT, "6");
    graph.setPathAttribute(7, GraphKmerAttribute.READ_COUNT, "7");
    graph.setPathAttribute(8, GraphKmerAttribute.READ_COUNT, "8");
    // 4 + 8
    assertEquals(12, GraphToPlot.overlappingReadCount(graph, 1, 2, 3));

    // 1 + 3 + 4 + 8
    assertEquals(16, GraphToPlot.overlappingReadCount(graph, 1, 2));

    // 3 + 4
    assertEquals(7, GraphToPlot.overlappingReadCount(graph, 4, 1, 2));

    // 7
    assertEquals(7, GraphToPlot.overlappingReadCount(graph, 1, 5));

    // 0
    assertEquals(0, GraphToPlot.overlappingReadCount(graph, 6, 1, 5));

    // 2 + 4 + 8
    assertEquals(14, GraphToPlot.overlappingReadCount(2, graph));

  }
  static final String EXPECTED_PALINDROME = ""
      + "digraph contigGraph {" + LS
      + "graph [rankdir=LR, ratio=fill]" + LS
      + "node0 [label=\"-2\\n5   3\"];" + LS
      + "node1 [label=\"+1\\n4   2\", color=red, shape=box];" + LS
      + "node2 [label=\"+3\\n6   4\"];" + LS
      + "node1 -> node0[color=red];" + LS
      + "node0 -> node1[color=red];" + LS
      + "node2 -> node0[];" + LS
      + "node0 -> node2[color=red];" + LS
      + "node2 -> node0[color=red];" + LS
      + "node1 -> node0[];" + LS
      + "}" + LS;

  public void testPalindrome() throws IOException {
    final MemoryPrintStream graphOut = new MemoryPrintStream();
    final MemoryPrintStream linksOut = new MemoryPrintStream();
    final Map<String, String> contigAttr = new HashMap<>();
    final Map<String, String> pathAttr = new HashMap<>();
    contigAttr.put(GraphKmerAttribute.READ_COUNT, "desc");
    pathAttr.put(GraphKmerAttribute.READ_COUNT, "desc");
    final GraphKmerAttribute graph = new GraphKmerAttribute(2, contigAttr, pathAttr);
    for (int i = 0; i < 3; ++i) {
      final int length = i + 4;
      contig(graph, length, length * (i + 1));
    }
    graph.addPath(new PathArray(1, 2));
    graph.addPath(new PathArray(1, -2));
    graph.addPath(new PathArray(2, 3));
    graph.addPath(new PathArray(2, -3));
    graph.setPathAttribute(1, GraphKmerAttribute.READ_COUNT, "2");
    graph.setContigAttribute(1, GraphKmerAttribute.READ_COUNT, "3");
    GraphToPlot.writeDotFile(graph, graphOut.printStream(), linksOut.printStream(), 2, 1L);
    mNano.check("graph-to-plot-palindrome.txt", graphOut.toString());
    //assertEquals(EXPECTED_PALINDROME, graphOut.toString());
  }
  static final String EXPECTED_CROSS = ""
      + "digraph contigGraph {" + LS
      + "graph [rankdir=LR, ratio=fill]" + LS
      + "node0 [label=\"+1\\n3   0\", color=red, shape=box];" + LS
      + "node1 [label=\"+2\\n3   0\"];" + LS
      + "node2 [label=\"+3\\n3   0\"];" + LS
      + "node3 [label=\"+4\\n3   0\"];" + LS
      + "node4 [label=\"+5\\n3   0\"];" + LS
      + "node1 -> node4[];" + LS
      + "node0 -> node1[label=\"2\"];" + LS
      + "node2 -> node1[];" + LS
      + "node1 -> node3[label=\"2\"];" + LS
      + "}" + LS;

  public void testCross() throws IOException {
    final MemoryPrintStream graphOut = new MemoryPrintStream();
    final MemoryPrintStream linksOut = new MemoryPrintStream();
    final Map<String, String> contigAttr = new HashMap<>();
    final Map<String, String> pathAttr = new HashMap<>();
    contigAttr.put(GraphKmerAttribute.READ_COUNT, "desc");
    pathAttr.put(GraphKmerAttribute.READ_COUNT, "desc");
    final MutableGraph graph = GraphMapCliTest.makeGraph(0, new String[]{"AAA", "AAA", "AAA", "AAA", "AAA"}, new long[][]{{1, 2, 4}, {3, 2, 5}}, contigAttr, pathAttr);
    graph.setPathAttribute(1, GraphKmerAttribute.READ_COUNT, "2");
    graph.setContigAttribute(1, GraphKmerAttribute.READ_COUNT, "3");
    GraphToPlot.writeDotFile(graph, graphOut.printStream(), linksOut.printStream(), 2, 1L);
    mNano.check("graph-to-plot-cross.txt", graphOut.toString());
    //assertEquals(EXPECTED_CROSS, graphOut.toString());
  }

  public static MutableGraph makeGraph(String[] contigs, long[][] paths) {
    final Map<String, String> attributes = new HashMap<>();
    attributes.put(GraphKmerAttribute.READ_COUNT, "read count");
    final MutableGraph graph = new GraphKmerAttribute(2, attributes, attributes);
    for (String c : contigs) {
      graph.addContig(new ContigString(c));
    }
    for (long[] path : paths) {
      graph.addPath(new PathArray(path));
    }
    return graph;
  }

  @Override
  protected AbstractCli getCli() {
    return new GraphToPlot();
  }
}
