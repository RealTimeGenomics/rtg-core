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

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Map;
import java.util.UUID;

import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.assembler.graph.io.GraphWriter;
import com.rtg.launcher.AbstractParamsCliTest;
import com.rtg.launcher.ParamsCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.store.StoreDirProxy;
import com.rtg.util.test.FileHelper;

/**
 */
public class GraphMapCliTest extends AbstractParamsCliTest<GraphMapParams> {
  public static GraphKmerAttribute makeGraph(int contigOverlap, String[] contigs, long[][] paths) {
    return makeGraph(contigOverlap, contigs, paths, Collections.<String, String>emptyMap(), Collections.<String, String>emptyMap());
  }
  public static GraphKmerAttribute makeGraph(int contigOverlap, String[] contigs, long[][] paths, Map<String, String> contigAttributes, Map<String, String> pathAttributes) {
    final GraphKmerAttribute graph = new GraphKmerAttribute(contigOverlap, contigAttributes, pathAttributes);
    graph.addContigAttribute(GraphKmerAttribute.READ_COUNT, "foo");
    graph.addPathAttribute(GraphKmerAttribute.READ_COUNT, "foo");
    for (String c : contigs) {
      graph.addContig(new ContigString(c));
    }
    for (long[] path : paths) {
      graph.addPath(new PathArray(path));
    }
    return graph;
  }
  public void testFlags() throws IOException {
    final MutableGraph builtGraph = makeGraph(29
        , new String[] {"ACGT", "GGGG", "TTAA"}
        , new long[][] {{1, 2}, {3, 2}});
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final File graph = new File(tmpDir, "input");
      final File reads = ReaderTestUtils.getDNADir(">a" + StringUtils.LS + "ACGTACGTACGTACGTACGTACGTACGTAC" + StringUtils.LS);
      assertTrue(graph.mkdir());
      try {
        GraphWriter.write(builtGraph, new StoreDirProxy(graph), "monkey", Collections.<UUID>emptySet());

        try {
          final File fileList = new File(tmpDir, "fileList");
          FileUtils.stringToFile(reads.toString() + StringUtils.LS, fileList);

          checkHandleFlagsErr();
          final File output = new File(tmpDir, "bar");
          String err = checkHandleFlagsErr("-o", output.toString(), reads.toString());
          TestUtils.containsAll(err, "You must provide a value for -g DIR");
          err = checkHandleFlagsErr("-o", output.toString(), "-g", graph.toString());
          TestUtils.containsAll(err, "No input files specified");
          err = checkHandleFlagsErr(reads.toString(), "-g", graph.toString());
          TestUtils.containsAll(err, "You must provide a value for -o DIR");

          err = checkHandleFlagsErr("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-w", "-1");
          TestUtils.containsAll(err, "--word", "at least 1");
          err = checkHandleFlagsErr("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-s", "-1");
          TestUtils.containsAll(err, "--step", "at least 1");

          err = checkHandleFlagsErr("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-a", "-1");
          TestUtils.containsAll(err, "The specified flag \"--mismatches\" has invalid value \"-1\"");

          err = checkHandleFlagsErr("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-m", "15", "-M", "10");
          TestUtils.containsAll(err, "--max-insert should be larger than --min-insert");

          checkHandleFlags("-o", output.toString(), reads.toString(), "-g", graph.toString(), "-a", "10%");
          checkHandleFlags("-o", output.toString(), reads.toString(), "-g", graph.toString());
          FileHelper.deleteAll(output);
          checkMainInitOk("-o", output.toString(), "-I", fileList.toString(), "-g", graph.toString());

          final File inputGraph =  new File(tmpDir, "graph");
          assertTrue(inputGraph.mkdir());
          GraphWriter.write(new GraphKmerAttribute(29), new StoreDirProxy(inputGraph), "foo" , Collections.<UUID>emptySet());
          final CFlags flags =  new CFlags("foo", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
          GraphMapCli.initLocalFlags(flags);
          flags.setFlags("-o", output.toString(), "-g", inputGraph.toString(), "-I", fileList.toString(), "-w", "29", "-s", "12");
          GraphMapParams params = GraphMapCli.makeParamsLocal(flags);
          assertEquals(29, params.wordSize());
          assertEquals(reads, params.reads().get(0));
          assertEquals(12, params.stepSize());
          flags.setFlags("-o", output.toString(), "-I", fileList.toString(), "-g", inputGraph.toString());
          params = GraphMapCli.makeParamsLocal(flags);
          assertEquals(18, params.wordSize());
          assertEquals(18, params.stepSize());
        } finally {
          FileHelper.deleteAll(graph);
        }

      } finally {
        FileHelper.deleteAll(reads);

      }
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testInitParams() {
    checkHelp("graphmap [OPTION]... -g DIR -o DIR DIR+",
        "-g DIR -o DIR -I FILE",
        "-g, --graph=DIR ", "graph of the assembly to map against",
        "-w, --word=INT", "word size",
        "-s, --step=INT", "step size",
        "-a, --mismatches=INT", "number of bases that may mismatch in an alignment",
        "-m, --min-insert=INT", "minimum insert size between fragments",
        "-M, --max-insert=INT", "maximum insert size between fragments",
        "-f, --454", "SDF containing 454 reads",
        "-F, --input-list-454", "file containing",
        "-j, --mate-pair", "SDF containing mate pair reads",
        "-J, --input-list-mate-pair"
    );
  }

  @Override
  protected ParamsCli<GraphMapParams> getParamsCli() {
    return new GraphMapCli();
  }
}
