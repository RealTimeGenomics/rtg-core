/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
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
