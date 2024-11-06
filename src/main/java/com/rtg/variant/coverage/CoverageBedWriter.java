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
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import com.rtg.bed.BedRecord;
import com.rtg.bed.BedWriter;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.Utils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.intervals.IntervalComparator;

/**
 * Handles BED and BEDGRAPH output formats
 */
public class CoverageBedWriter extends CoverageProcessor {

  private static final int COVERAGE_DP = GlobalFlags.getIntegerValue(CoreGlobalFlags.COVERAGE_DP);

  static final String DEFAULT_LABEL = "coverage";

  private final BedWriter mOut;
  private final CoverageParams mParams;
  private final boolean mSort;
  private final List<BedRecord> mBuffered = new ArrayList<>();
  private String mLastSequenceName = null;
  private String mBedLabel = DEFAULT_LABEL;

  /**
   * Create a new coverage writer
   * @param params output options
   * @throws IOException if there is a problem creating the output stream
   */
  public CoverageBedWriter(CoverageParams params) throws IOException {
    this(params.bedStream(), params);
  }

  /**
   * Create a new coverage writer
   * @param out send output here
   * @param params output options
   */
  public CoverageBedWriter(OutputStream out, CoverageParams params) {
    mOut = new BedWriter(out);
    mSort = params.perRegion(); // Per-region can result in out-of-order records, so we need to sort
    mParams = params;
  }

  @Override
  public void init() throws IOException {
    if (mParams.bedgraphOutput()) {
      mOut.writeln(VERSION_STRING + ", Coverage BEDGRAPH output " + COVERAGE_OUTPUT_VERSION);
    } else {
      mOut.writeln(VERSION_STRING + ", Coverage BED output " + COVERAGE_OUTPUT_VERSION);
    }
    if (CommandLine.getCommandLine() != null) {
      mOut.writeComment("CL" + TB + CommandLine.getCommandLine());
    }
    mOut.writeComment("RUN-ID" + TB + CommandLine.getRunId());
    if (mParams.bedgraphOutput()) {
      mOut.writeln("#sequence\tstart\tend\tcoverage");
      mOut.writeln("track type=bedGraph name=coverage");
    } else {
      mOut.writeln("#sequence\tstart\tend\tlabel\tcoverage");
    }
  }

  @Override
  public void finalCoverageRegion(String name, int start, int end, double coverage) throws IOException {
    final BedRecord rec;
    if (mParams.bedgraphOutput()) {
      rec = new BedRecord(name, start, end, Utils.realFormat(coverage, COVERAGE_DP));
    } else {
      rec = new BedRecord(name, start, end, mBedLabel, Utils.realFormat(coverage, COVERAGE_DP));
    }
    if (mSort) {
      if (!name.equals(mLastSequenceName)) {
        flush();
        mLastSequenceName = name;
      }
      mBuffered.add(rec);
    } else {
      mOut.write(rec);
    }
  }

  private void flush() throws IOException {
    mBuffered.sort(IntervalComparator.SINGLETON);
    for (BedRecord rec : mBuffered) {
      mOut.write(rec);
    }
    mBuffered.clear();
  }

  @Override
  public void setRegionLabel(String label) {
    mBedLabel = label.trim();
    if (mBedLabel.length() == 0) {
      mBedLabel = DEFAULT_LABEL;
    }
  }

  @Override
  public void close() throws IOException {
    flush();
    mOut.close();
  }
}
