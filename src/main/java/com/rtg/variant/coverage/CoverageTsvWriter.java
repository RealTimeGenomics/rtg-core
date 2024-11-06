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

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.Environment;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;

/**
 * Handles per-base output format.
 */
public class CoverageTsvWriter extends CoverageProcessor {

  /** Coverage file version string */
  public static final String VERSION_STRING = "#Version " + Environment.getVersion();

  static final String COVERAGE_OUTPUT_VERSION = "v1.1";
  private static final String TB = "\t";

  private final OutputStream mOut;

  /**
   * Create a new coverage writer
   * @param params output options
   * @throws IOException if there is a problem creating the output stream
   */
  public CoverageTsvWriter(CoverageParams params) throws IOException {
    this(params.bedStream());
  }

  /**
   * Create a new coverage writer
   * @param out send output here
   */
  public CoverageTsvWriter(OutputStream out) {
    mOut = out;
  }

  @Override
  public void init() throws IOException {
    mOut.write((VERSION_STRING + ", Coverage output " + COVERAGE_OUTPUT_VERSION + LS).getBytes());
    if (CommandLine.getCommandLine() != null) {
      mOut.write(("#CL" + TB + CommandLine.getCommandLine() + LS).getBytes());
    }
    mOut.write(("#RUN-ID" + TB + CommandLine.getRunId() + LS).getBytes());
    mOut.write(("#sequence\tposition\tunique-count\tambiguous-count\tscore" + LS).getBytes());
  }

  @Override
  public void finalCoveragePosition(String name, int position, int ih1, int ihgt1, double coverage) throws IOException {
    mOut.write((name + TB + position + TB + ih1 + TB + ihgt1 + TB + Utils.realFormat(coverage, 2) + LS).getBytes());
  }

  @Override
  public void close() throws IOException {
    mOut.close();
  }
}
