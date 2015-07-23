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

    int mMinimumCoverageForBreadth = 1;
    int mSmoothing = 0;
    boolean mErrorRates = false;
    boolean mTsvOutput = false;
    boolean mBedgraphOutput = false;
    boolean mOnlyMappedRegions = false;
    boolean mOutputIndex = true;
    int mChunkSize = 10000;
    boolean mDisableHtmlReport = false;

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
     * Only produce output for regions which have coverage, but still tracks stats for those regions
     * @param onlyMapped true to only produce output for regions with coverage.
     * @return this builder, so calls can be chained.
     */
    public CoverageParamsBuilder onlyMappedRegions(final boolean onlyMapped) {
      mOnlyMappedRegions = onlyMapped;
      return self();
    }

    /**
     * Minimum coverage for breadth calculation
     *
     * @param minCoverage minimum coverage
     * @return this builder, so calls can be chained.
     */
    public  CoverageParamsBuilder minimumCoverageForBreadth(int minCoverage) {
      mMinimumCoverageForBreadth = minCoverage;
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
  private final boolean mTsvOutput;
  private final boolean mBedgraphOutput;
  private final boolean mOnlyMappedRegions;
  private final int mMinimumCoverageForBreadth;
  private final int mChunkSize;
  private final boolean mDisableHtmlReport;

  private final boolean mOutputIndex;

  /**
   * @param builder the builder object.
   */
  CoverageParams(final CoverageParamsBuilder builder) {
    super(builder);
    mSmoothing = builder.mSmoothing;
    mErrorRates = builder.mErrorRates;
    mTsvOutput = builder.mTsvOutput;
    mBedgraphOutput = builder.mBedgraphOutput && !mTsvOutput; // TSV takes priority if both set
    mOnlyMappedRegions = builder.mOnlyMappedRegions;
    mMinimumCoverageForBreadth = builder.mMinimumCoverageForBreadth;
    mOutputIndex = builder.mOutputIndex;
    mChunkSize = builder.mChunkSize;
    mDisableHtmlReport = builder.mDisableHtmlReport;
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
   * Get whether to only report stats on mapped regions.
   * @return true means only write out statistics for regions with coverage.
   */
  public boolean onlyMappedRegions() {
    return mOnlyMappedRegions;
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
   * @return minimum coverage for breadth calculations
   */
  public int minimumCoverageForBreadth() {
    return mMinimumCoverageForBreadth;
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

  @Override
  public String toString() {
    return "CoverageParams"
      + " mapped reads=" + mapped().toString()
      + " smoothing=" + mSmoothing
      + " error rates=" + Boolean.valueOf(mErrorRates).toString()
      + " minCoverageForBreadth=" + mMinimumCoverageForBreadth
      + LS + super.toString();
  }
}
