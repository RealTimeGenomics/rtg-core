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
import java.util.Collections;

import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.GraphKmerAttribute;
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
public class PacBioCliTest extends AbstractParamsCliTest<PacBioParams> {
  @Override
  protected ParamsCli<PacBioParams> getParamsCli() {
    return new PacBioCli();
  }

  public void testFlags() throws IOException {
    final MutableGraph builtGraph = GraphMapCliTest.makeGraph(29
        , new String[] {"ACGT", "GGGG", "TTAA"}
        , new long[][] {{1, 2}, {3, 2}});
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final File graph = new File(tmpDir, "input");
      final File reads = ReaderTestUtils.getDNADir(">a" + StringUtils.LS + "ACGTACGTACGTACGTACGTACGTACGTAC" + StringUtils.LS);
      assertTrue(graph.mkdir());
      try {
        GraphWriter.write(builtGraph, new StoreDirProxy(graph), "monkey", Collections.emptySet());

        try {
          final File fileList = new File(tmpDir, "fileList");
          FileUtils.stringToFile(reads.toString() + StringUtils.LS, fileList);

          checkHandleFlagsErr();
          final File output = new File(tmpDir, "bar");
          String err = checkHandleFlagsErr("-o", output.toString(), reads.toString());
          assertTrue(err, err.contains("You must provide a value for -g DIR"));
          err = checkHandleFlagsErr("-o", output.toString(), "-g", graph.toString());
          assertTrue(err, err.contains("No input files specified"));
          err = checkHandleFlagsErr(reads.toString(), "-g", graph.toString());
          assertTrue(err, err.contains("You must provide a value for -o DIR"));

          checkHandleFlags("-o", output.toString(), reads.toString(), "-g", graph.toString());
          FileHelper.deleteAll(output);
          checkMainInitOk("-o", output.toString(), "-I", fileList.toString(), "-g", graph.toString());

          final File inputGraph =  new File(tmpDir, "graph");
          assertTrue(inputGraph.mkdir());
          GraphWriter.write(new GraphKmerAttribute(29), new StoreDirProxy(inputGraph), "foo" , Collections.emptySet());
          final CFlags flags =  new CFlags("foo", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
          PacBioCli.initLocalFlags(flags);
          flags.setFlags("-o", output.toString(), "-g", inputGraph.toString(), "-I", fileList.toString());
          final PacBioParams params = PacBioCli.makeParamsLocal(flags);
          assertEquals(reads, params.reads().get(0));
          assertEquals(output, params.directory());
          assertEquals(inputGraph, params.graph());
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

  public void testHelp() {
    checkHelp("addpacbio [OPTION]... -g DIR -o DIR SDF+",
        "-I, --input-list-file=FILE ", "file containing a list of SDF directories",
        "-g, --graph=DIR ", "graph of the assembly to map against",
        "SDF directories containing reads to map"
    );
  }
}
