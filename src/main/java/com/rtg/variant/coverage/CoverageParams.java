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

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.sam.SingleMappedParams;
import com.rtg.util.test.params.ParamsNoField;

/**
 */
public final class CoverageParams extends SingleMappedParams {

  /** name of bed graph file */
  public static final String BEDGRAPH_NAME = "coverage.bedgraph";

  /** name of bed file */
  public static final String BED_NAME = "coverage.bed";

  /** name of tab separated file */
  public static final String TSV_NAME = "coverage.tsv";

  /**
   * Creates a <code>CoverageParamsBuilder</code>.
   * @return the builder.
   */
  public static CoverageParamsBuilder builder() {
    return new CoverageParamsBuilder();
  }

  /**
   * A builder class for <code>CoverageParams</code>.
   */
  public static final class CoverageParamsBuilder extends SingleMappedParamsBuilder<CoverageParamsBuilder> {

    int mFoldTargetPercent = 80;
    int mMinimumCoverageThreshold = 1;
    int mSmoothing = 0;
    boolean mErrorRates = false;
    boolean mTsvOutput = false;
    boolean mBedgraphOutput = false;
    private boolean mPerRegion = false;
    boolean mOutputIndex = true;
    int mChunkSize = 10000;
    boolean mDisableHtmlReport = false;
    private boolean mBinarizeBed = false;
    private boolean mIncludeDeletions = false;

    @Override
    protected CoverageParamsBuilder self() {
      return this;
    }

    /**
     * Sets the smoothing level.
     *
     * @param neighbours the number of neighbours on each side to become smooth with.
     * @return this builder, so calls can be chained.
     */
    public CoverageParamsBuilder smoothing(final int neighbours) {
      mSmoothing = neighbours;
      return self();
    }

    /**
     * Turns on reporting of sequencer error rates.
     *
     * @param errorRates true means generate a report.
     * @return this builder, so calls can be chained.
     */
    public CoverageParamsBuilder errorRates(final boolean errorRates) {
      mErrorRates = errorRates;
      return self();
    }

    /**
     * Turns on inclusion of deletion operations in coverage values
     * @param includeDeletions true if deletions should be counted
     * @return this builder, so calls can be chained.
     */
    public CoverageParamsBuilder includeDeletions(final boolean includeDeletions) {
      mIncludeDeletions = includeDeletions;
      return self();
    }

    /**
     * Turns on coverage binarization
     * @param binarize true if BED output should be binarized
     * @return this builder, so calls can be chained.
     */
    public CoverageParamsBuilder binarizeBed(final boolean binarize) {
      mBinarizeBed = binarize;
      return self();
    }

    /**
     * Turns on per-region aggregation
     * @param perRegion true if BED output should be aggregate per region
     * @return this builder, so calls can be chained.
     */
    public CoverageParamsBuilder perRegion(final boolean perRegion) {
      mPerRegion = perRegion;
      return self();
    }

    /**
     * Turns on tab separated value output.
     * @param tsv true means generate output in tab separated value format.
     * @return this builder, so calls can be chained.
     */
    public CoverageParamsBuilder tsvOutput(final boolean tsv) {
      mTsvOutput = tsv;
      return self();
    }

    /**
     * Turns on BEDGRAPH output.
     * @param bedgraph true means generate output in the BEDGRAPH format.
     * @return this builder, so calls can be chained.
     */
    public CoverageParamsBuilder bedgraphOutput(final boolean bedgraph) {
      mBedgraphOutput = bedgraph;
      return self();
    }

    /**
     * Minimum coverage for breadth calculation
     *
     * @param minCoverage minimum coverage
     * @return this builder, so calls can be chained.
     */
    public  CoverageParamsBuilder minimumCoverageThreshold(int minCoverage) {
      mMinimumCoverageThreshold = minCoverage;
      return self();
    }

    /**
     * Target percentage coverage for fold penalty statistic
     *
     * @param foldPct target percentage coverage
     * @return this builder, so calls can be chained.
     */
    public CoverageParamsBuilder foldTargetPercent(int foldPct) {
      mFoldTargetPercent = foldPct;
      return self();
    }

    /**
     * @param val true if TABIX index should be created for output coverage files.
     * @return this builder, so calls can be chained.
     */
    public CoverageParamsBuilder outputIndex(boolean val) {
      mOutputIndex = val;
      return self();
    }

    /**
     * Disable the output of HTML reports
     * @param value true to disable the HTML reports
     * @return this builder
     */
    public CoverageParamsBuilder disableHtmlReport(boolean value) {
      mDisableHtmlReport = value;
      return self();
    }

    /**
     * Creates a CoverageParams using the current builder
     * configuration.
     * @return the new CoverageParams
     */
    public CoverageParams create() {
      return new CoverageParams(this);
    }
  }

  private final int mSmoothing;
  private final boolean mErrorRates;
  private final boolean mBinarizeBed;
  private final boolean mIncludeDeletions;
  private final boolean mPerRegion;
  private final boolean mTsvOutput;
  private final boolean mBedgraphOutput;
  private final int mMinimumCoverageThreshold;
  private final int mChunkSize;
  private final boolean mDisableHtmlReport;
  private final int mFoldTargetPercent;
  private final boolean mOutputIndex;

  /**
   * @param builder the builder object.
   */
  CoverageParams(final CoverageParamsBuilder builder) {
    super(builder);
    mSmoothing = builder.mSmoothing;
    mErrorRates = builder.mErrorRates;
    mBinarizeBed = builder.mBinarizeBed;
    mIncludeDeletions = builder.mIncludeDeletions;
    mTsvOutput = builder.mTsvOutput;
    mPerRegion = builder.mPerRegion && !mTsvOutput;
    mBedgraphOutput = builder.mBedgraphOutput && !mTsvOutput; // TSV takes priority if both set
    mMinimumCoverageThreshold = builder.mMinimumCoverageThreshold;
    mOutputIndex = builder.mOutputIndex;
    mChunkSize = builder.mChunkSize;
    mDisableHtmlReport = builder.mDisableHtmlReport;
    mFoldTargetPercent = builder.mFoldTargetPercent;
  }

  /**
   *  Get smoothing level.
   * @return smoothing level.
   */
  public int smoothing() {
    return mSmoothing;
  }

  /**
   *  Get error rates flag.
   * @return true means generate a sequencer error rates report.
   */
  public boolean errorRates() {
    return mErrorRates;
  }

  /**
   * Get the output file name
   * @return the output file name
   */
  @ParamsNoField
  public String fileName() {
    return mTsvOutput ? TSV_NAME : mBedgraphOutput ? BEDGRAPH_NAME : BED_NAME;
  }

  /**
   * @return the stream for writing the coverage bed file.
   * @throws IOException whenever.
   */
  @ParamsNoField
  public OutputStream bedStream() throws IOException {
    return outStream(fileName());
  }

  /**
   * @return output file to be used
   */
  @ParamsNoField
  public File outFile() {
    return outFile(fileName());
  }

  /**
   * @return true if deletion operators should be included in coverage counts (e.g. for variant callability)
   */
  public boolean includeDeletions() {
    return mIncludeDeletions;
  }

  /**
   * @return true if the BED should be binarized
   */
  public boolean binarizeBed() {
    return mBinarizeBed;
  }

  /**
   * @return true if BED/BEDGRAPH should report mean per region
   */
  public boolean perRegion() {
    return mPerRegion;
  }

  /**
   * Get whether tab separated value is being output.
   * @return true means generate the tab separated value format.
   */
  public boolean tsvOutput() {
    return mTsvOutput;
  }

  /**
   * @return true if output is being written as <code>BED</code> format
   */
  @ParamsNoField
  public boolean bedOutput() {
    return !mBedgraphOutput && !mTsvOutput;
  }

  /**
   * @return true if output is being written as <code>BEDGRAPH</code> format
   */
  public boolean bedgraphOutput() {
    return mBedgraphOutput;
  }

  /**
   * @return minimum coverage for breadth calculations and binarization
   */
  public int minimumCoverageThreshold() {
    return mMinimumCoverageThreshold;
  }

  /**
   * @return true if TABIX index should be output.
   */
  public boolean outputIndex() {
    return mOutputIndex;
  }

  /**
   * @return the number of positions to include when dividing work into chunks.
   */
  public int chunkSize() {
    return mChunkSize;
  }

  /** @return true if HTML report output should be disabled */
  public boolean disableHtmlReport() {
    return mDisableHtmlReport;
  }

  /**
   * @return the target percentage of bases when computing fold penalty statistic.
   */
  public int foldTargetPercent() {
    return mFoldTargetPercent; 
  }

  @Override
  public String toString() {
    return "CoverageParams"
      + " mapped reads=" + mapped()
      + " smoothing=" + mSmoothing
      + " error rates=" + mErrorRates
      + " min coverage threshold=" + mMinimumCoverageThreshold
      + " binarize bed=" + mBinarizeBed
      + " include deletions=" + mIncludeDeletions
      + LS + super.toString();
  }
}
