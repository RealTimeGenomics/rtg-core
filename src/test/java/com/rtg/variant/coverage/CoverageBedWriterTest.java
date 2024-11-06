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

package com.rtg.variant.coverage;

import java.io.IOException;

import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

public class CoverageBedWriterTest extends TestCase {

  private static final String GRAPH_LINE = "track type=bedGraph name=coverage";

  public void testBedGraphWriter() throws IOException {
    final MemoryPrintStream out = new MemoryPrintStream();
    final CoverageParams params = CoverageParams.builder().bedgraphOutput(true).create();
    final CoverageBedWriter cw = new CoverageBedWriter(out.outputStream(), params);
    cw.init();
    cw.finalCoverageRegion("hello", 1, 20, 3);
    TestUtils.containsAll(out.toString(),
        CoverageBedWriter.VERSION_STRING,
        "#RUN-ID",
        "Coverage BEDGRAPH output " + CoverageBedWriter.COVERAGE_OUTPUT_VERSION,
        GRAPH_LINE,
        "hello\t1\t20\t3"
    );
  }

  public void testBedWriter() throws IOException {
    final String[] args = CommandLine.getCommandLine() != null ? CommandLine.getCommandArgs() : null;
    try {
    CommandLine.setCommandArgs("foo", "bar");
    final MemoryPrintStream out = new MemoryPrintStream();
    final CoverageParams params = CoverageParams.builder().create();
    final CoverageBedWriter cw = new CoverageBedWriter(out.outputStream(), params);
    cw.init();
    cw.finalCoverageRegion("hello", 1, 20, 3);
    final String actual = out.toString();
    TestUtils.containsAll(actual,
        CoverageBedWriter.VERSION_STRING,
        "#RUN-ID",
        "#CL\tfoo bar",
        "Coverage BED output " + CoverageBedWriter.COVERAGE_OUTPUT_VERSION,
        "hello\t1\t20\tcoverage\t3"
    );
    assertFalse(actual.contains(GRAPH_LINE));
    } finally {
      if (args != null) {
        CommandLine.setCommandArgs(args);
      }
    }
  }
}
