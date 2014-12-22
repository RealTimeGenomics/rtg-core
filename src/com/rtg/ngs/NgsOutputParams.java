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
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

import net.sf.samtools.SAMReadGroupRecord;

/**
 * Holds all the parameters needed for doing a build and a search.
 *
 */
public class NgsOutputParams extends IntegralAbstract {

  /** repeat file name */
  public static final String REPEATS_FILE_NAME = "ambiguous";

  /** unmapped file name */
  public static final String UNMAPPED_FILE_NAME = "unmapped";

  /** SAM single-end alignments file name */
  public static final String ALIGNMENTS_SAM_FILE_NAME = "alignments.sam";

  /** BAM single-end alignments file name */
  public static final String ALIGNMENTS_BAM_FILE_NAME = "alignments.bam";

  /** SAM mated alignments file name */
  public static final String MATED_SAM_FILE_NAME = "mated.sam";

  /** BAM mated alignments file name */
  public static final String MATED_BAM_FILE_NAME = "mated.bam";

  /** SAM unmated alignments file name */
  public static final String UNMATED_SAM_FILE_NAME = "unmated.sam";

  /** BAM unmated alignments file name */
  public static final String UNMATED_BAM_FILE_NAME = "unmated.bam";

  /** SAM unmapped file name */
  public static final String UNMAPPED_SAM_FILE_NAME = "unmapped.sam";

  /** BAM unmapped file name */
  public static final String UNMAPPED_BAM_FILE_NAME = "unmapped.bam";

  /**
   * Creates a NgsOutputParams builder.
   * @return the builder.
   */
  public static NgsOutputParamsBuilder builder() {
    return new NgsOutputParamsBuilder();
  }

  private final boolean mProgress;

  private final File mOutputDir;

  private final File mTempFilesDir;

  private final int mNumOutputFiles;

  private final boolean mTabular;

  private final boolean mSorted;

  private final boolean mBam;

  private final boolean mSam;

  private final boolean mSdf;

  private final boolean mUnify;

  private final NgsFilterParams mFilterParams;

  private final NgsResultStreamHandler mResultStreamHandler;

  private final boolean mKeepIntermediate;

  private final boolean mMergeMatchResults;

  private final boolean mMergeAlignmentResults;

  private final boolean mOutputUnmated;

  private final boolean mOutputUnmapped;

  private final boolean mOutputReadNames;

  private final boolean mOutputProteinSequences;

  private final SAMReadGroupRecord mReadGroupRecord;

  private final boolean mCalibrate;

  private final boolean mSvPrep;

  private final boolean mOutputIndex;

  private final ReferenceRegions mCalibrateRegions;

  /**
   * @param builder the builder object.
   */
  public NgsOutputParams(final NgsOutputParamsBuilder builder) {
    mProgress = builder.mProgress;
    mOutputDir = builder.mOutputDir;
    mTempFilesDir = builder.mTempFilesDir;
    mNumOutputFiles = builder.mNumOutputFiles;
    mTabular = builder.mTabular;
    mSorted = builder.mSorted;
    mBam = builder.mBam;
    mSam = builder.mSam;
    mSdf = builder.mSdf;
    mUnify = builder.mUnify;
    mKeepIntermediate = builder.mKeepIntermediate;
    mMergeMatchResults = builder.mMergeMatchResults;
    mMergeAlignmentResults = builder.mMergeAlignmentResults;
    mOutputUnmated = builder.mOutputUnmated;
    mOutputUnmapped = builder.mOutputUnmapped;
    mOutputReadNames = builder.mOutputReadNames;
    mOutputProteinSequences = builder.mOutputProteinSequences;
    mFilterParams = builder.mFilterParams;
    mReadGroupRecord = builder.mReadGroupRecord;
    mResultStreamHandler = new NgsResultStreamHandler(mOutputDir, mTempFilesDir, mFilterParams.zip());
    mCalibrate = builder.mCalibrate;
    mCalibrateRegions = builder.mCalibrateRegions;
    mSvPrep = builder.mSvPrep;
    mOutputIndex = builder.mOutputIndex;
  }

  /**
   * Get the progress flag.
   *
   * @return the progress flag. (true iff progress is to be output).
   */
  public boolean progress() {
    return mProgress;
  }

  /**
   * Get the filter parameters.
   * @return the filter parameters.
   */
  public NgsFilterParams filter() {
    return mFilterParams;
  }

  /**
   * Returns the current output filter to be used
   * @return output filter
   */
  public OutputFilter outFilter() {
    return mFilterParams.outputFilter();
  }

  /**
   * Get the name of a child file in the output directory where all results are placed.
   * @param child the name of the child.
   * @return the name of the file.
   */
  public File file(final String child) {
    return mResultStreamHandler.file(child);
  }

  /**
   * Get a stream to the output file.
   * @return the stream.
   * @throws IOException if the stream cannot be created
   */
  public OutputStream outStream() throws IOException {
    return mResultStreamHandler.outStream();
  }

  /**
   * Get a stream to the output file.
   * @return the stream.
   * @throws IOException if the stream cannot be created
   */
  public OutputStream matedSamStream() throws IOException {
    return mResultStreamHandler.matedSamStream();
  }

  /**
   * Return the result stream handler
   * @return as said
   */
  public NgsResultStreamHandler resultStreamHandler() {
    return mResultStreamHandler;
  }

  /**
   * Stream to write repeats to.
   * @return the stream
   * @throws IOException if the stream cannot be created
   */
  public OutputStream repeatsStream() throws IOException {
    return mResultStreamHandler.repeatStream();
  }

  /**
   * Stream to write unmapped reads to.
   * @return the stream
   * @throws IOException if the stream cannot be created
   */
  public OutputStream unmappedStream() throws IOException {
    return mResultStreamHandler.unmappedStream();
  }

  /**
   * Get the error limit. Used in output to limit what is written out.
   *
   * @return the error limit (&gte;0).
   */
  public int errorLimit() {
    return mFilterParams.errorLimit();
  }

  /**
   * Get the number of repeated hits to be remembered for each read.
   * @return the number of repeated hits to be remembered for each read.
   */
  public int topN() {
    return mFilterParams.topN();
  }

  /**
   * Get the number of repeated hits to be remembered for each read.
   * @return the number of repeated hits to be remembered for each read.
   */
  public int maxTopResults() {
    return mFilterParams.maxTopResults();
  }

  /**
   * Get the output directory.
   * @return the output directory.
   */
  public File directory() {
      return mOutputDir;
  }

  /**
   * Get the directory for temporary files.
   * @return if the temporary files directory is null, the normal output directory is returned.
   */
  public File tempFilesDirectory() {
    if (mTempFilesDir != null) {
      return mTempFilesDir;
    } else {
      return mOutputDir;
    }
  }

  /**
   * Get the number of output files to produce.
   * @return number of output files
   */
  public int numberOfOutputFiles() {
    return mNumOutputFiles;
  }

  /**
   * Returns whether compressing output files is enabled.
   *
   * @return compress output
   */
  public boolean isCompressOutput() {
    return filter().zip();
  }

  /**
   * Check if reads with repeated hits should be excluded.
   * @return true if reads with repeated hits should be excluded.
   */
  public boolean exclude() {
    return mFilterParams.exclude();
  }

  /**
   * Should unmated results be output in the paired end case
   * @return the truth
   */
  public boolean outputUnmated() {
    return mOutputUnmated;
  }

  /**
   * Should unmapped results be output
   * @return the truth
   */
  public boolean outputUnmapped() {
    return mOutputUnmapped;
  }

  /**
   * Check if numeric identifiers should be used to name sequences (rather than the names of the sequences
   * as specified in the input files).
   * @return true iff numeric identifiers should be used.
   */
  public boolean useids() {
    return mFilterParams.useids();
  }

  @Override
  public int hashCode() {
    return Utils.hash(new Object[] {mProgress, mOutputDir, mTempFilesDir, mFilterParams, mTabular, mSorted, mBam, mSam, mSdf, mKeepIntermediate, mMergeAlignmentResults, mMergeMatchResults, mNumOutputFiles});
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final NgsOutputParams that = (NgsOutputParams) obj;
    return Utils.equals(
        new Object[] {this.mProgress, this.mFilterParams, this.mOutputDir, this.mTempFilesDir, this.mTabular, this.mSorted, this.mBam, this.mSam, this.mSdf, this.mNumOutputFiles, this.mKeepIntermediate, this.mMergeMatchResults, this.mMergeAlignmentResults},
        new Object[] {that.mProgress, that.mFilterParams, that.mOutputDir, that.mTempFilesDir, that.mTabular, that.mSorted, that.mBam, that.mSam, that.mSdf, that.mNumOutputFiles, that.mKeepIntermediate, that.mMergeMatchResults, that.mMergeAlignmentResults}
    );
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("" + "NgsOutputParams" + " progress=").append(mProgress).append(" output=").append(mOutputDir).append(" filterParams={").append(mFilterParams).append("}").append(" tabular=").append(mTabular).append(" sorted=").append(mSorted).append(" tempDir=").append(mTempFilesDir).append(" numoutfiles=").append(mNumOutputFiles).append(" bam=").append(mBam).append(" sam=").append(mSam).append(" sdf=").append(mSdf).append(" keepIntermediate=").append(mKeepIntermediate).append(" mergeMatchResults=").append(mMergeMatchResults).append(" mergeAlignmentResults=").append(mMergeAlignmentResults);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mOutputDir != null /*&& mOutputDir.isDirectory()*/); //extra check isnt a good idea given our current testing
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    Exam.globalIntegrity(mFilterParams);
    return true;
  }

  /**
   * @return If the output should be tabular
   */
  public boolean tabular() {
    return mTabular;
  }

  /**
   * @return if output to be sorted
   */
  public boolean sorted() {
    return mSorted;
  }

  /**
   * whether output should be BAM
   * @return true if output should be BAM
   */
  public boolean bam() {
    return mBam;
  }

  /**
   * whether should output SAM
   * @return true if should output SAM
   */
  public boolean sam() {
    return mSam;
  }

  /**
   * whether should output SDF
   * @return true if should output SDF
   */
  public boolean sdf() {
    return mSdf;
  }

  /**
   * whether should unify SAM/BAM output
   * @return true if should
   */
  public boolean unify() {
    return mUnify;
  }

  /**
   * @return whether intermediate files should be kept
   */
  public boolean keepIntermediate() {
    return mKeepIntermediate;
  }

  /**
   * @return whether intermediate match output should be merged
   */
  public boolean mergeMatchResults() {
    return mMergeMatchResults;
  }

  /**
   * @return whether intermediate alignment output should be merged
   */
  public boolean mergeAlignmentResults() {
    return mMergeAlignmentResults;
  }

  /**
   * @return whether to use read names in output instead of read ids.
   */
  public boolean outputReadNames() {
    return mOutputReadNames;
  }

  /**
   * @return whether to output protein sequence information in <code>mapx</code> output.
   */
  public boolean outputProteinSequences() {
    return mOutputProteinSequences;
  }

  /**
   * @return {@link SAMReadGroupRecord} for current run, can be null
   */
  public SAMReadGroupRecord readGroup() {
    return mReadGroupRecord;
  }

  /**
   * @return true if calibration files need to be produced, false otherwise
   */
  public boolean calibrate() {
    return mCalibrate;
  }

  /**
   * @return the bed regions to use if calibrating on a subset of the whole genome
   */
  public ReferenceRegions calibrateRegions() {
    return mCalibrateRegions;
  }

  /**
   * @return true if structural variation preparation set is to be performed.
   */
  public boolean svprep() {
    return mSvPrep;
  }

  /**
   * @return true if TABIX or BAM index should be output.
   */
  public boolean outputIndex() {
    return mOutputIndex;
  }

}
