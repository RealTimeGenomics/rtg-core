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
package com.rtg.variant.eval;

import java.io.File;

import com.rtg.launcher.OutputModuleParams;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.vcf.VcfUtils;

/**
 * Params object for VcfEvalTask
 */
public final class VcfEvalParams extends OutputModuleParams {

  /**
   * Creates a <code>VcfEvalParamsBuilder</code> builder.
   * @return the builder.
   */
  public static VcfEvalParamsBuilder builder() {
    return new VcfEvalParamsBuilder();
  }

  /**
   * A builder class for <code>VcfEvalParams</code>.
   */
  public static final class VcfEvalParamsBuilder extends OutputModuleParamsBuilder<VcfEvalParamsBuilder> {
    File mBaselineFile;
    File mCallsFile;
    File mTemplateFile;
    RocSortOrder mSortOrder = RocSortOrder.DESCENDING;
    String mScoreField = VcfUtils.FORMAT_GENOTYPE_QUALITY;
    String mSampleName;
    boolean mUseAllRecords = false;
    int mNumberThreads = 1;
    private boolean mSquashPloidy = false;
    int mMaxLength = -1;
    boolean mRtgStats = false;
    boolean mOutputBaselineTp = false;
    boolean mOutputSlopeFiles = false;
    private RegionRestriction mRestriction = null;
    private File mBedRegionsFile = null;

    @Override
    protected VcfEvalParamsBuilder self() {
      return this;
    }

    /**
     * Sets the baseline file.
     *
     * @param baseline the baseline file.
     * @return this builder, so calls can be chained.
     */
    public VcfEvalParamsBuilder baseLineFile(final File baseline) {
      mBaselineFile = baseline;
      return self();
    }

    /**
     * Sets the calls file.
     *
     * @param calls the calls file.
     * @return this builder, so calls can be chained.
     */
    public VcfEvalParamsBuilder callsFile(final File calls) {
      mCallsFile = calls;
      return self();
    }

    /**
     * Sets the template file.
     *
     * @param template the template file.
     * @return this builder, so calls can be chained.
     */
    public VcfEvalParamsBuilder templateFile(final File template) {
      mTemplateFile = template;
      return self();
    }

    /**
     * Set to use all records from VCF, by default we only use PASS records.
     *
     * @param useAllRecords set to <code>true</code> to use all records for ROC.
     * @return this builder, so calls can be chained.
     */
    public VcfEvalParamsBuilder useAllRecords(final boolean useAllRecords) {
      mUseAllRecords = useAllRecords;
      return self();
    }

    /**
     * Sets the sample name to use in multiple sample VCF input.
     *
     * @param sampleName the sample name to use in multiple sample VCF input.
     * @return this builder, so calls can be chained.
     */
    public VcfEvalParamsBuilder sampleName(final String sampleName) {
      mSampleName = sampleName;
      return self();
    }

    /**
     * Sets the VCF ROC scoring field.
     *
     * @param scoreField the VCF format field to use as the ROC score.
     * @return this builder, so calls can be chained.
     */
    public VcfEvalParamsBuilder scoreField(final String scoreField) {
      mScoreField = scoreField;
      return self();
    }

    /**
     * Sets the sort order for the ROC score.
     *
     * @param sortOrder the sort order for the ROC score.
     * @return this builder, so calls can be chained.
     */
    public VcfEvalParamsBuilder sortOrder(final RocSortOrder sortOrder) {
      mSortOrder = sortOrder;
      return self();
    }

    /**
     * @param value number of threads to use when evaluating variants
     * @return this builder, so calls can be chained
     */
    public VcfEvalParamsBuilder numberThreads(int value) {
      mNumberThreads = value;
      return self();
    }

    /**
     * @param squashPloidy true if heterozygous diploid calls should be squashed to homozygous alternative
     * @return this builder, so calls can be chained
     */
    public VcfEvalParamsBuilder squashPloidy(boolean squashPloidy) {
      mSquashPloidy = squashPloidy;
      return self();
    }

    /**
     * @param maxLength the maximum length variant to consider
     * @return this builder, so calls can be chained
     */
    public VcfEvalParamsBuilder maxLength(int maxLength) {
      mMaxLength = maxLength;
      return self();
    }


    /**
     * @param rtgStats should RTG specific output be produced
     * @return this builder, so calls can be chained
     */
    public VcfEvalParamsBuilder rtgStats(boolean rtgStats) {
      mRtgStats = rtgStats;
      return self();
    }

    /**
     * @param outputSlopeFiles if set, output the files for ROC slope analysis
     * @return this builder, so calls can be chained
     */
    public VcfEvalParamsBuilder outputSlopeFiles(boolean outputSlopeFiles) {
      mOutputSlopeFiles = outputSlopeFiles;
      return self();
    }

    /**
     * @param outputBaselineTp if set, output the baseline version of the true positives file too
     * @return this builder, so calls can be chained
     */
    public VcfEvalParamsBuilder outputBaselineTp(boolean outputBaselineTp) {
      mOutputBaselineTp = outputBaselineTp;
      return self();
    }

    /**
     * Sets a restriction on the records that will be processed. The format is a reference sequence name,
     * followed by an optional range specification. Only records that match the reference name and start
     * within the range (if specified) will be processed. A name of null indicates no filtering.
     * @param restriction a reference sequence name
     * @return this builder, so calls can be chained.
     */
    public VcfEvalParamsBuilder restriction(final RegionRestriction restriction) {
      mRestriction = restriction;
      return this;
    }

    /**
     * Set the bed file to use which specifies regions
     * @param bedRegionsFile the bed file which specifies regions
     * @return this builder
     */
    public VcfEvalParamsBuilder bedRegionsFile(File bedRegionsFile) {
      mBedRegionsFile = bedRegionsFile;
      return this;
    }

    /**
     * Creates a <code>VcfEvalParams</code> using the current builder
     * configuration.
     * @return the new <code>VcfEvalParams</code>
     */
    public VcfEvalParams create() {
      return new VcfEvalParams(this);
    }
  }

  private final File mBaselineFile;
  private final File mCallsFile;
  private final File mTemplateFile;
  private final RegionRestriction mRestriction;
  private final File mBedRegionsFile;
  private final String mScoreField;
  private final RocSortOrder mSortOrder;
  private final String mSampleName;
  private final int mNumberThreads;
  private final boolean mUseAllRecords;
  private final boolean mSquashPloidy;
  private final int mMaxLength;
  private final boolean mRtgStats;
  private final boolean mOutputBaselineTp;
  private final boolean mOutputSlopeFiles;


  /**
   * @param builder the builder object.
   */
  protected VcfEvalParams(VcfEvalParamsBuilder builder) {
    super(builder);
    mBaselineFile = builder.mBaselineFile;
    mCallsFile = builder.mCallsFile;
    mTemplateFile = builder.mTemplateFile;
    mRestriction = builder.mRestriction;
    mBedRegionsFile = builder.mBedRegionsFile;
    mSortOrder = builder.mSortOrder;
    mScoreField = builder.mScoreField;
    mSampleName = builder.mSampleName;
    mNumberThreads = builder.mNumberThreads;
    mUseAllRecords = builder.mUseAllRecords;
    mSquashPloidy = builder.mSquashPloidy;
    mMaxLength = builder.mMaxLength;
    mRtgStats = builder.mRtgStats;
    mOutputSlopeFiles = builder.mOutputSlopeFiles;
    mOutputBaselineTp = builder.mOutputBaselineTp;
  }

  /**
   * Get the baseline file.
   * @return the baseline file.
   */
  public File baselineFile() {
    return mBaselineFile;
  }

  /**
   * Get the calls file.
   * @return the calls file.
   */
  public File callsFile() {
    return mCallsFile;
  }

  /**
   * Get the template file.
   * @return the template file.
   */
  public File templateFile() {
    return mTemplateFile;
  }


  /**
   * @return The region to restrict iteration of the VCF records to.
   */
  public RegionRestriction restriction() {
    return mRestriction;
  }

  /**
   * @return a bed file containing the regions to process, or null for no bed region based filtering.
   */
  public File bedRegionsFile() {
    return mBedRegionsFile;
  }

  /**
   * Get the sort field from VCF to use for ROC curves.
   * @return the sort field from VCF format field to use for ROC curves.
   */
  public String scoreField() {
    return mScoreField;
  }

  /**
   * Get the sort order to use for ROC curves.
   * @return the sort order to use in ROC curves.
   */
  public RocSortOrder sortOrder() {
    return mSortOrder;
  }

  /**
   * Get the sample name for multiple sample VCF input.
   * @return the sample name to use for multiple sample VCF input.
   */
  public String sampleName() {
    return mSampleName;
  }

  /**
   * Get whether to use all VCF records for ROC or not
   * @return <code>true</code> if true all VCF records are used, if <code>false</code>, only PASS VCF records will be used.
   */
  public boolean useAllRecords() {
    return mUseAllRecords;
  }

  /**
   * @return true if heterozygous diploid calls should be squashed to homozygous alternative
   */
  public boolean squashPloidy() {
    return mSquashPloidy;
  }

  /**
   * @return the number of threads to run in
   */
  public int numberThreads() {
    return mNumberThreads;
  }

  /**
   * @return the limit on variant length
   */
  public int maxLength() {
    return mMaxLength;
  }
  /**
   * @return true if RTG specific stats should be included in output
   */
  public boolean rtgStats() {
    return mRtgStats;
  }

  /**
   * @return true if slope files should be written
   */
  public boolean outputSlopeFiles() {
    return mOutputSlopeFiles;
  }



  /**
   * @return true if the baseline version of the true positives file should be output
   */
  public boolean outputBaselineTp() {
    return mOutputBaselineTp;
  }

  @Override
  public String toString() {
    return "Baseline file=" + mBaselineFile.getPath() + ", Calls file=" + mCallsFile.getPath() + ", Template file=" + mTemplateFile.getPath() + ", score field=" + mScoreField + ", sort order=" + mSortOrder + ", sample name=" + mSampleName + ", num threads=" + mNumberThreads + ", use all records=" + mUseAllRecords + ", squash ploidy=" + mSquashPloidy + ", max length=" + mMaxLength + ", rtg stats=" + mRtgStats + ", baseline tp=" + mOutputBaselineTp + ", output params=" + super.toString();
  }

}
