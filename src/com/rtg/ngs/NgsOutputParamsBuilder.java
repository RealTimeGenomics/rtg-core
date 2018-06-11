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
package com.rtg.ngs;

import java.io.File;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.intervals.ReferenceRegions;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 * A builder class for <code>NgsOutputParams</code>.
 */
@TestClass("com.rtg.ngs.NgsOutputParamsTest")
public class NgsOutputParamsBuilder {

  protected boolean mProgress = false;

  protected File mOutputDir = null;

  protected File mTempFilesDir = null;

  protected int mNumOutputFiles = 1;

  protected boolean mTabular = false;

  protected boolean mSorted = false;

  protected boolean mBam = false;

  protected boolean mUnify = false;

  protected boolean mSam = true;

  protected boolean mSdf = false;

  protected boolean mMergeMatchResults = true;

  protected boolean mMergeAlignmentResults = true;

  protected boolean mOutputUnmated = false;

  protected boolean mOutputUnmapped = false;

  protected boolean mCalibrate = false;

  protected boolean mSvPrep = false;

  protected NgsFilterParams mFilterParams = NgsFilterParams.builder()
    .outputFilter(OutputFilter.TOPN_PAIRED_END)
    .useids(true)
    .errorLimit(MapFlags.MAX_SCORE)
    .create();

  protected boolean mIgnoreShort = false;

  //this is to set if we dont want to delete intermediate files
  protected boolean mKeepIntermediate = false;

  protected boolean mOutputReadNames = false;

  protected boolean mOutputProteinSequences = true;

  protected boolean mOutputIndex = true;

  //store read group record to be added in final SAM file
  protected SAMReadGroupRecord mReadGroupRecord = null;

  protected ReferenceRegions mCalibrateRegions = null;

  /**
   * @param progress turn progress on or off
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder progress(final boolean progress) {
    mProgress = progress;
    return this;
  }

  /**
   * @param dir output directory
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder outputDir(final File dir) {
    mOutputDir = dir;
    return this;
  }

  /**
   * @param dir directory for temporary files (may be null)
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder tempFilesDir(final File dir) {
    mTempFilesDir = dir;
    return this;
  }

  /**
   * @param n number of alignment output files to produce
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder numberOfOutputFiles(final int n) {
    mNumOutputFiles = n;
    return this;
  }

  /**
   * @param tabular tabular output
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder tabular(final boolean tabular) {
    mTabular = tabular;
    return this;
  }

  /**
   * @param sorted sorted output
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder sorted(final boolean sorted) {
    mSorted = sorted;
    return this;
  }

  /**
   * @param params filter params
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder filterParams(final NgsFilterParams params) {
    mFilterParams = params;
    return this;
  }

  /**
   * Write alignments in BAM
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder bam(final boolean val) {
    mBam = val;
    return this;
  }

  /**
   * Write alignments in SAM
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder sam(final boolean val) {
    mSam = val;
    return this;
  }

  /**
   * Write alignments in SDF
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder sdf(final boolean val) {
    mSdf = val;
    return this;
  }

  /**
   * Unify SAM/BAM output
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder unify(final boolean val) {
    mUnify = val;
    return this;
  }

  /**
   * Ignore reads shorter than index word size from output and statistics
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder ignoreShort(final boolean val) {
    mIgnoreShort = val;
    return this;
  }

  /**
   * Do not delete intermediate files
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder keepIntermediateFiles(final boolean val) {
    mKeepIntermediate = val;
    return this;
  }

  /**
   * merge results of match phase, default = true
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder mergeMatchResults(final boolean val) {
    mMergeMatchResults = val;
    return this;
  }

  /**
   * merge results of alignment phase, default = true
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder mergeAlignmentResults(final boolean val) {
    mMergeAlignmentResults = val;
    return this;
  }

  /**
   * Output unmated results, default = false
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder outputUnmated(final boolean val) {
    mOutputUnmated = val;
    return this;
  }

  /**
   * Output unmapped results, default = false
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder outputUnmapped(final boolean val) {
    mOutputUnmapped = val;
    return this;
  }

  /**
   * Output read names, default = false
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder outputReadNames(boolean val) {
    mOutputReadNames  = val;
    return this;
  }

  /**
   * Output protein sequences, default = true
   * @param val the value
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder outputProteinSequences(boolean val) {
    mOutputProteinSequences  = val;
    return this;
  }

  /**
   * add read group to output file
   * @param samReadGroupRecord SAM read group record
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder readGroup(SAMReadGroupRecord samReadGroupRecord) {
    mReadGroupRecord = samReadGroupRecord;
    return this;
  }

  /**
   * set if calibration files need to be produced (default is true)
   * @param val true if calibration needs to be performed, false otherwise
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder calibrate(boolean val) {
    mCalibrate = val;
    return this;
  }

  /**
   * set if you'd like calibration to be restricted to the regions provided in the bed file
   * @param val a bed file
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder calibrateRegions(ReferenceRegions val) {
    mCalibrateRegions = val;
    return this;
  }

  /**
   * set if svprep stage is to be run
   * @param val true if svprep is to be performed, false otherwise
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder svprep(boolean val) {
    mSvPrep = val;
    return this;
  }

  /**
   * @param val true if TABIX or BAM index should be created for output SAM or BAM files.
   * @return this builder, so calls can be chained.
   */
  public NgsOutputParamsBuilder outputIndex(boolean val) {
    mOutputIndex = val;
    return this;
  }

  /**
   * Creates a NgsParams using the current builder
   * configuration.
   * @return the new NgsParams
   */
  public NgsOutputParams create() {
    return new NgsOutputParams(this);
  }
}
