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

package com.rtg.variant.cnv.segment;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.OutputModuleParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.util.test.params.ParamsNoField;

/** Parameters for CNV segmentation */
public class SegmentParams extends OutputModuleParams {

  /** @return a builder instance. */
  public static SegmentParamsBuilder builder() {
    return new SegmentParamsBuilder();
  }

  /** The parameter builder class. */
  public static class SegmentParamsBuilder extends OutputModuleParamsBuilder<SegmentParamsBuilder> {

    private ReaderParams mGenome;
    private File mSummaryRegionsFile;
    private File mPanelFile;
    private File mControlFile;
    private File mCaseFile;
    private int mPrecomputedColumn = -1;
    private String mSampleName = "SAMPLE";
    private String mCoverageColumnName = "coverage";
    private String mPanelCoverageColumnName = "normalized-coverage";
    private int mGcBins = 10;
    private int mMinBins = 1;
    private int mMinSegments = 1;
    private int mMaxSegments = Integer.MAX_VALUE;
    private Double mMinNormControlCoverage = 0.1;
    private Double mMinControlCoverage;
    private Double mMinCaseCoverage;
    private double mAleph;
    private double mAlpha;
    private double mBeta;
    private double mMinLogR = 0.2;
    private boolean mAbsorbSingletons;
    private boolean mGraphviz;

    /**
     * Sets the genome reader.
     * @param params the parameters used for output.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder genome(final ReaderParams params) {
      mGenome = params;
      return self();
    }

    /**
     * Sets a BED file containing regions for reporting CNV interactions.
     * @param file the BED file.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder summaryRegionsFile(File file) {
      mSummaryRegionsFile = file;
      return this;
    }

    /**
     * Sets the BED file supplying panel of normals coverage data.
     * @param file the BED file.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder panelFile(File file) {
      mPanelFile = file;
      return this;
    }

    /**
     * Sets the BED file supplying control sample coverage data.
     * @param file the BED file.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder controlFile(File file) {
      mControlFile = file;
      return this;
    }

    /**
     * Sets the BED file supplying case sample coverage data.
     * @param file the BED file.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder caseFile(File file) {
      mCaseFile = file;
      return this;
    }

    /**
     * Sets the sample name for output
     * @param name the sample name
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder sampleName(String name) {
      mSampleName = name;
      return this;
    }

    /**
     * Sets the index of an input column to segment directly. The default is to
     * perform normalization and ratio computation from input coverage columns.
     * @param precomputedColumn direct segmentation column index.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder precomputedColumn(int precomputedColumn) {
      mPrecomputedColumn = precomputedColumn;
      return this;
    }

    /**
     * Sets the name of the input data column containing coverage statistics.
     * @param name the column name
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder coverageColumnName(String name) {
      mCoverageColumnName = name;
      return this;
    }

    /**
     * Sets the name of the input panel of normals data column containing normalized coverage statistics.
     * @param name the column name
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder panelCoverageColumnName(String name) {
      mPanelCoverageColumnName = name;
      return this;
    }

    /**
     * Sets the number of bins used for GC normalization.
     * @param gcBins the number of bins
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder gcBins(int gcBins) {
      mGcBins = gcBins;
      return this;
    }

    /**
     * Sets the minimum number of bins a segment may contain
     * @param bins the number of bins
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder minBins(int bins) {
      mMinBins = bins;
      return this;
    }

    /**
     * Sets the minimum number of segments a chromosome may be reduced to.
     * @param segments the number of segments.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder minSegments(int segments) {
      mMinSegments = segments;
      return this;
    }

    /**
     * Sets the maximum number of segments a chromosome may be reduced to.
     * @param segments the number of segments.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder maxSegments(int segments) {
      mMaxSegments = segments;
      return this;
    }

    /**
     * Sets the minimum normalized control/panel coverage for a region to be included.
     * @param value the minimum normalized control coverage.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder minNormControlCoverage(Double value) {
      mMinNormControlCoverage = value;
      return this;
    }

    /**
     * Sets the minimum control sample coverage for a region to be included.
     * @param coverage the minimum coverage.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder minControlCoverage(Double coverage) {
      mMinControlCoverage = coverage;
      return this;
    }

    /**
     * Sets the minimum case sample coverage for a region to be included.
     * @param coverage the minimum coverage.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder minCaseCoverage(Double coverage) {
      mMinCaseCoverage = coverage;
      return this;
    }

    /**
     * Sets the inter-segment distance weighting factor
     * @param aleph the weight.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder aleph(double aleph) {
      mAleph = aleph;
      return this;
    }

    /**
     * Sets the intra-segment distance weighting factor
     * @param alpha the weight.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder alpha(double alpha) {
      mAlpha = alpha;
      return this;
    }

    /**
     * Sets the segmentation sensitivity factor. Higher values produce fewer segments.
     * @param beta the segmentation sensitivity.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder beta(double beta) {
      mBeta = beta;
      return this;
    }

    /**
     * Sets minimum absolute log R to call a deviation as a CNV.
     * @param logR the minimum absolute log R.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder minLogR(double logR) {
      mMinLogR = logR;
      return this;
    }

    /**
     * If set, absorb single-bin segments into the closest neighbour.
     * @param value true if singletons should be absorbed.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder absorbSingletons(boolean value) {
      mAbsorbSingletons = value;
      return this;
    }

    /**
     * If set, produce Graphviz visualization output.
     * @param value true to output Graphviz visualization.
     * @return this builder, so calls can be chained.
     */
    public SegmentParamsBuilder graphviz(boolean value) {
      mGraphviz = value;
      return this;
    }

    @Override
    protected SegmentParamsBuilder self() {
      return this;
    }

    /** @return a new parameter object from current settings */
    public SegmentParams create() {
      return new SegmentParams(this);
    }
  }


  private final ReaderParams mGenome;
  private final File mPanelFile;
  private final File mControlFile;
  private final File mCaseFile;
  private final File mSummaryRegionsFile;
  private final String mSampleName;
  private final String mCoverageColumnName;
  private final String mPanelCoverageColumnName;
  private final int mPrecomputedColumn;
  private final int mGcBins;
  private final int mMinBins;
  private final int mMinSegments;
  private final int mMaxSegments;
  private final double mMinLogR;
  private final Double mMinNormControlCoverage;
  private final Double mMinControlCoverage;
  private final Double mMinCaseCoverage;
  private final double mAleph;
  private final double mAlpha;
  private final double mBeta;
  private final boolean mAbsorbSingletons;
  private final boolean mGraphviz;

  /** @param builder the builder object. */
  protected SegmentParams(SegmentParamsBuilder builder) {
    super(builder);
    mGenome = builder.mGenome;
    mPanelFile = builder.mPanelFile;
    mControlFile = builder.mControlFile;
    mCaseFile = builder.mCaseFile;
    mSummaryRegionsFile = builder.mSummaryRegionsFile;
    mSampleName = builder.mSampleName;
    mPrecomputedColumn = builder.mPrecomputedColumn;
    mCoverageColumnName = builder.mCoverageColumnName;
    mPanelCoverageColumnName = builder.mPanelCoverageColumnName;
    mGcBins = builder.mGcBins;
    mMinBins = builder.mMinBins;
    mMinNormControlCoverage = builder.mMinNormControlCoverage;
    mMinControlCoverage = builder.mMinControlCoverage;
    mMinCaseCoverage = builder.mMinCaseCoverage;
    mAleph = builder.mAleph;
    mAlpha = builder.mAlpha;
    mBeta = builder.mBeta;
    mMinLogR = builder.mMinLogR;
    mMinSegments = builder.mMinSegments;
    mMaxSegments = builder.mMaxSegments;
    mAbsorbSingletons = builder.mAbsorbSingletons;
    mGraphviz = builder.mGraphviz;
  }

  /** @return parameters for the reference genome. */
  public ReaderParams genome() {
    return mGenome;
  }

  /** @return the BED file supplying panel of normals coverage data, or null. */
  public File panelFile() {
    return mPanelFile;
  }

  /** @return the BED file supplying control sample coverage data, or null. */
  public File controlFile() {
    return mControlFile;
  }

  /** @return the BED file supplying case sample coverage data, or null. */
  public File caseFile() {
    return mCaseFile;
  }

  /** @return the sample name for output. */
  public String sampleName() {
    return mSampleName;
  }

  /** @return a BED file containing regions for reporting CNV interactions, or null. */
  public File summaryRegionsFile() {
    return mSummaryRegionsFile;
  }

  /** @return the index of an input column to segment directly, or -1. */
  public int precomputedColumn() {
    return mPrecomputedColumn;
  }

  /** @return the name of the input data column containing coverage statistics. */
  public String coverageColumnName() {
    return mCoverageColumnName;
  }

  /** @return the name of the input panel of normals data column containing normalized coverage statistics. */
  public String panelCoverageColumnName() {
    return mPanelCoverageColumnName;
  }

  /** @return the number of bins used for GC normalization. */
  public int gcBins() {
    return mGcBins;
  }

  /** @return the minimum normalized control/panel coverage for a region to be included. */
  public Double minNormControlCoverage() {
    return mMinNormControlCoverage;
  }

  /** @return the minimum control sample coverage for a region to be included. */
  public Double minControlCoverage() {
    return mMinControlCoverage;
  }

  /** @return the minimum case sample coverage for a region to be included. */
  public Double minCaseCoverage() {
    return mMinCaseCoverage;
  }

  /** @return the inter-segment distance weighting factor */
  public double aleph() {
    return mAleph;
  }

  /** @return the intra-segment distance weighting factor */
  public double alpha() {
    return mAlpha;
  }

  /** @return the segmentation sensitivity factor. Higher values produce fewer segments. */
  public double beta() {
    return mBeta;
  }

  /** @return the minimum number of bins a segment may contain. */
  public int minBins() {
    return mMinBins;
  }

  /** @return the minimum absolute log R to call a deviation as a CNV. */
  public double minLogR() {
    return mMinLogR;
  }

  /** @return the minimum number of segments a chromosome may be reduced to. */
  public int minSegments() {
    return mMinSegments;
  }

  /** @return the maximum number of segments a chromosome may be reduced to. */
  public int maxSegments() {
    return mMaxSegments;
  }

  /** @return true if singletons should be absorbed. */
  public boolean absorbSingletons() {
    return mAbsorbSingletons;
  }

  /** @return true if Graphviz output should be produced. */
  public boolean graphviz() {
    return mGraphviz;
  }

  @Override
  @ParamsNoField
  public void close() throws IOException {
    if (mGenome != null) {
      mGenome.close();
    }
  }
}

