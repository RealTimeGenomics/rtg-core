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

package com.rtg.variant.coverage;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.Environment;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;

/**
 *         Date: 16/02/12
 *         Time: 10:01 AM
 */
public class CoverageWriter implements CoverageProcessor {
  private static final String TB = "\t";
  /** Coverage file version string */
  public static final String VERSION_STRING = "#Version " + Environment.getVersion();
  static final String DEFAULT_LABEL = "coverage";

  private final OutputStream mOut;
  private final CoverageParams mParams;
  private static final String COVERAGE_OUTPUT_VERSION = "v1.0";

  private String mBedLabel = DEFAULT_LABEL;

  /**
   * Create a new coverage writer
   * @param out send output here
   * @param params output options
   */
  public CoverageWriter(OutputStream out, CoverageParams params) {
    mOut = out;
    mParams = params;
  }

  @Override
  public void init() throws IOException {
    if (mParams.bedOutput()) {
      mOut.write((VERSION_STRING + ", Coverage BED output " + COVERAGE_OUTPUT_VERSION + LS).getBytes());
    } else if (mParams.bedgraphOutput()) {
      mOut.write((VERSION_STRING + ", Coverage BEDGRAPH output " + COVERAGE_OUTPUT_VERSION + LS).getBytes());
    } else {
      mOut.write((VERSION_STRING + ", Coverage output " + COVERAGE_OUTPUT_VERSION + LS)
          .getBytes());
    }
    if (CommandLine.getCommandLine() != null) {
      mOut.write(("#CL" + TB + CommandLine.getCommandLine() + LS).getBytes());
    }
    //mOut.write()
    mOut.write(("#RUN-ID" + TB + CommandLine.getRunId() + LS).getBytes());
    if (mParams.tsvOutput()) {
      mOut.write(("#sequence\tposition\tunique-count\tambiguous-count\tscore" + LS).getBytes());
    } else if (mParams.bedgraphOutput()) {
      mOut.write(("#sequence\tstart\tend\tcoverage" + LS).getBytes());
      mOut.write(("track type=bedGraph name=coverage" + LS).getBytes());
    } else {
      mOut.write(("#sequence\tstart\tend\tlabel\tcoverage" + LS).getBytes());
    }
  }

  @Override
  public void finalCoveragePosition(String name, int position, int ih1, int ihgt1, double coverage) throws IOException {
    mOut.write((name + TB + position + TB + ih1 + TB + ihgt1 + TB + Utils.realFormat(coverage, 2) + LS).getBytes());
  }

  @Override
  public void finalCoverageRegion(String name, int start, int end, int coverage) throws IOException {
    if (mParams.bedgraphOutput()) {
      mOut.write((name + TB + start + TB + end + TB + coverage + LS).getBytes());
    } else {
      mOut.write((name + TB + start + TB + end + TB + mBedLabel + TB + coverage + LS).getBytes());
    }
  }

  public void setBedLabel(String label) {
    mBedLabel = label.trim();
    if (mBedLabel.length() == 0) {
      mBedLabel = DEFAULT_LABEL;
    }
  }

  @Override
  public void close() throws IOException {
    mOut.close();
  }
}
