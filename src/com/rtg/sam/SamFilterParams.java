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
package com.rtg.sam;

import java.io.File;

import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.InvalidParamsException;

/**
 * Parameters for selection of SAM records.
 */
public class SamFilterParams {

  /**
   * Creates a builder.
   * @return the builder.
   */
  public static SamFilterParamsBuilder builder() {
    return new SamFilterParamsBuilder();
  }

  /** The builder. */
  public static class SamFilterParamsBuilder {

    protected int mMaxAlignmentCount = -1;
    protected int mMinMapQ = -1;
    protected IntegerOrPercentage mMaxASMatedValue = null;
    protected IntegerOrPercentage mMaxASUnmatedValue = null;
    protected boolean mExcludeMated = false;
    protected boolean mExcludeUnmated = false;
    protected boolean mExcludeUnmapped = false;
    protected boolean mExcludeDuplicates = false;
    protected boolean mExcludeUnplaced = false;
    protected int mRequireUnsetFlags = 0;
    protected int mRequireSetFlags = 0;
    protected boolean mFindAndRemoveDuplicates = false; //this is for the detect and remove, rather than looking at the sam flag

    protected SamRegionRestriction mRestriction = null;
    protected File mBedRegionsFile = null;

    /**
     * SAM records having MAPQ less than this value will be filtered.
     * @param val the minimum MAPQ permitted or -1 for no filtering
     * @return this builder, so calls can be chained
     */
    public SamFilterParamsBuilder minMapQ(final int val) {
      mMinMapQ = val;
      return this;
    }

    /**
     * SAM records having alignment count over this value will be filtered.
     * @param val the maximum count or -1 for no filtering
     * @return this builder, so calls can be chained
     */
    public SamFilterParamsBuilder maxAlignmentCount(final int val) {
      mMaxAlignmentCount = val;
      return this;
    }

    /**
     * Mated SAM records having alignment score over this value will be filtered
     * @param val the maximum score as integer or percentage, or null for no filtering
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder maxMatedAlignmentScore(final IntegerOrPercentage val) {
      mMaxASMatedValue = val;
      return this;
    }

    /**
     * Unmated SAM records having alignment score over this value will be filtered
     * @param val the maximum score as integer or percentage, or null for no filtering
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder maxUnmatedAlignmentScore(final IntegerOrPercentage val) {
      mMaxASUnmatedValue = val;
      return this;
    }

    /**
     * Should unmated records be excluded. Default is false.
     * @param val true to exclude unmated records
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder excludeUnmated(final boolean val) {
      mExcludeUnmated = val;
      return this;
    }

    /**
     * Should mated records be excluded. Default is false.
     * @param val true to exclude mated records
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder excludeMated(final boolean val) {
      mExcludeMated = val;
      return this;
    }

    /**
     * Should unmapped records be excluded. Default is false.
     * @param val true to exclude unmapped records
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder excludeUnmapped(final boolean val) {
      mExcludeUnmapped = val;
      return this;
    }

    /**
     * Should records without an alignment start position be excluded. Default is false.
     * @param val true to exclude unmapped records
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder excludeUnplaced(final boolean val) {
      mExcludeUnplaced = val;
      return this;
    }

    /**
     * Should PCR and optical duplicate records be excluded. Default is false.
     * @param val true to exclude unmapped records
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder excludeDuplicates(final boolean val) {
      mExcludeDuplicates = val;
      return this;
    }

    /**
     * Specify mask indicating SAM flags that must be unset. Any record with any
     * of these flags set will be excluded. Default is not checking any of these flags.
     * @param flags mask indicating flags that must be unset.
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder requireUnsetFlags(final int flags) {
      mRequireUnsetFlags = flags;
      return this;
    }

    /**
     * Specify mask indicating SAM flags that must be set. Any record with any
     * of these flags unset will be excluded. Default is not checking any of these flags.
     * @param flags mask indicating flags that must be set.
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder requireSetFlags(final int flags) {
      mRequireSetFlags = flags;
      return this;
    }

    /**
     * Should duplicates be detected and removed. Default is true.
     * @param val true to exclude detected duplicates
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder findAndRemoveDuplicates(final boolean val) {
      mFindAndRemoveDuplicates = val;
      return this;
    }

    /**
     * Sets a restriction on the records that will be processed. The format is a reference sequence name,
     * followed by an optional range specification. Only records that match the reference name and start
     * within the range (if specified) will be processed. A name of null indicates no filtering.
     * @param restriction a reference sequence name
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder restriction(final String restriction) {
      if (restriction == null) {
        mRestriction = null;
      } else  {
        mRestriction = new SamRegionRestriction(restriction);
      }
      return this;
    }

    /**
     * Sets a restriction on the records that will be processed. The format is a reference sequence name,
     * followed by an optional range specification. Only records that match the reference name and start
     * within the range (if specified) will be processed. A name of null indicates no filtering.
     * @param restriction a reference sequence name
     * @return this builder, so calls can be chained.
     */
    public SamFilterParamsBuilder restriction(final SamRegionRestriction restriction) {
      mRestriction = restriction;
      return this;
    }

    /**
     * Set the bed file to use which specifies regions
     * @param bedRegionsFile the bed file which specifies regions
     * @return this builder
     */
    public SamFilterParamsBuilder bedRegionsFile(File bedRegionsFile) {
      mBedRegionsFile = bedRegionsFile;
      return this;
    }

    /**
     * Create the immutable version.
     * @return a <code>SamFilterParams</code> value
     */
    public SamFilterParams create() {
      return new SamFilterParams(this);
    }
  }

  private final int mMaxAlignmentCount;
  private final int mMinMapQ;
  private final IntegerOrPercentage mMaxASMatedValue;
  private final IntegerOrPercentage mMaxASUnmatedValue;
  private final int mRequireUnsetFlags;
  private final int mRequireSetFlags;
  private final boolean mExcludeUnmated;
  private final boolean mExcludeUnplaced;
  private final boolean mFindAndRemoveDuplicates; //detected version

  private final SamRegionRestriction mRestriction;
  private final File mBedRegionsFile;

  /**
   * @param builder the builder object.
   */
  public SamFilterParams(SamFilterParamsBuilder builder) {
    mMinMapQ = builder.mMinMapQ;
    mMaxAlignmentCount = builder.mMaxAlignmentCount;
    mMaxASMatedValue = builder.mMaxASMatedValue;
    mMaxASUnmatedValue = builder.mMaxASUnmatedValue;
    mExcludeUnplaced = builder.mExcludeUnplaced;

    int requireUnsetFlags = builder.mRequireUnsetFlags;
    if (builder.mExcludeDuplicates) {
      requireUnsetFlags |= SamBamConstants.SAM_PCR_OR_OPTICAL_DUPLICATE;
    }
    if (builder.mExcludeMated) {
      requireUnsetFlags |= SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR;
    }
    if (builder.mExcludeUnmapped) {
      requireUnsetFlags |= SamBamConstants.SAM_READ_IS_UNMAPPED;
    }
    mRequireUnsetFlags = requireUnsetFlags;
    mRequireSetFlags = builder.mRequireSetFlags;

    final int badFlags = mRequireUnsetFlags & mRequireSetFlags;
    if (badFlags != 0) {
      throw new InvalidParamsException("Conflicting SAM FLAG criteria that no record can meet: " + badFlags);
    }

    // This one is tricky - needs to hit reads that are !unmapped and !properpair
    mExcludeUnmated = builder.mExcludeUnmated;

    mRestriction = builder.mRestriction;
    mBedRegionsFile = builder.mBedRegionsFile;

    mFindAndRemoveDuplicates = builder.mFindAndRemoveDuplicates;
  }

  /**
   * Mated SAM records having alignment score over this value will be filtered.
   * @return the maximum score or null for no filtering
   */
  public IntegerOrPercentage maxMatedAlignmentScore() {
    return mMaxASMatedValue;
  }

  /**
   * Unmated SAM records having alignment score over this value will be filtered.
   * @return the maximum score or null for no filtering
   */
  public IntegerOrPercentage maxUnmatedAlignmentScore() {
    return mMaxASUnmatedValue;
  }

  /**
   * SAM records having MAPQ less than this value will be filtered.
   * @return the minimum MAPQ permitted or -1 for no filtering
   */
  public int minMapQ() {
    return mMinMapQ;
  }

  /**
   * SAM records having alignment count over this value will be filtered.
   * @return the maximum count or -1 for no filtering
   */
  public int maxAlignmentCount() {
    return mMaxAlignmentCount;
  }

  /**
   * True if unmated SAM records should be excluded.
   * @return exclusion status
   */
  public boolean excludeUnmated() {
    return mExcludeUnmated;
  }

  /**
   * True if SAM records without an alignment start position should be excluded.
   * @return exclusion status
   */
  public boolean excludeUnplaced() {
    return mExcludeUnplaced;
  }

  /**
   * A mask indicating SAM flags that must be unset. Any record with any
   * of these flags set will be excluded. Default is not checking any of these flags.
   * @return mask indicating flags that must be unset.
   */
  public int requireUnsetFlags() {
    return mRequireUnsetFlags;
  }

  /**
   * A mask indicating SAM flags that must be set. Any record with any
   * of these flags unset will be excluded. Default is not checking any of these flags.
   * @return mask indicating flags that must be set.
   */
  public int requireSetFlags() {
    return mRequireSetFlags;
  }

  /**
   * Should duplicates be detected and removed.
   * @return true to exclude detected duplicates
   */
  public boolean findAndRemoveDuplicates() {
    return mFindAndRemoveDuplicates;
  }

  /**
   * @return The region to restrict iteration of the SAM record to.
   */
  public SamRegionRestriction restriction() {
    return mRestriction;
  }

  /**
   * Returns reference sequence name to filter on.  Null means no filtering.
   * @return reference sequence name to use
   */
  public String restrictionTemplate() {
    return mRestriction == null ? null : mRestriction.getSequenceName();
  }

  /**
   * Gets the zero based start position of restricted calling. -1 means no position filtering.
   *
   * @return the zero based restriction start position inclusive.
   */
  public int restrictionStart() {
    return mRestriction == null ? -1 : mRestriction.getStart();
  }

  /**
   * Gets the zero based exclusive end position of restricted calling.
   *
   * @return the zero based exclusive restriction end position.
   */
  public int restrictionEnd() {
    return mRestriction == null ? -1 : mRestriction.getEnd();
  }


  /**
   * @return a bed file containing the regions to process, or null for no bed region based filtering.
   */
  public File bedRegionsFile() {
    return mBedRegionsFile;
  }

  @Override
  public String toString() {
    return "SamFilterParams"
    + " minMapQ=" + minMapQ()
    + " maxAlignmentCount=" + maxAlignmentCount()
    + " maxMatedAlignmentScore=" + maxMatedAlignmentScore()
    + " maxUnmatedAlignmentScore=" + maxUnmatedAlignmentScore()
    + " excludeUnmated=" + excludeUnmated()
    + " excludeUnplaced=" + excludeUnplaced()
    + " requireSetFlags=" + requireSetFlags()
    + " requireUnsetFlags=" + requireUnsetFlags()
    + " regionTemplate=" + restrictionTemplate();
  }

}
