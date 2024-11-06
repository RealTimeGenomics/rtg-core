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
package com.rtg.ngs;

import java.io.Serializable;

import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public final class NgsFilterParams extends IntegralAbstract implements Serializable {

  /** maximum alignment score allowed for the mated output */
  public static final String MAX_MATED_MISMATCH_THRESHOLD = "10%";

  /** maximum alignment score allowed for the unmated output */
  public static final String MAX_UNMATED_MISMATCH_THRESHOLD = "10%";
  /**
   * Creates a NgsFilterParams builder.
   * @return the builder.
   */
  public static NgsFilterParamsBuilder builder() {
    return new NgsFilterParamsBuilder();
  }

  /**
   * A builder class for <code>NgsFilterParams</code>.
   */
  public static class NgsFilterParamsBuilder {

    protected OutputFilter mOutputFilter = OutputFilter.TOPN_PAIRED_END;

    protected boolean mZip = false;

    protected int mTopN = 5;

    protected int mMaxTopResults = 5;

    protected boolean mUseSequenceIds = false;

    protected boolean mExcludeRepeats = false;

    protected int mErrorLimit = 3;

    protected IntegerOrPercentage mUnmatedMaxMismatches = IntegerOrPercentage.valueOf(MAX_UNMATED_MISMATCH_THRESHOLD);
    protected IntegerOrPercentage mMatedMaxMismatches = IntegerOrPercentage.valueOf(MAX_MATED_MISMATCH_THRESHOLD);

    protected double mMaxEScore = 10.0;
    protected double mMinBitScore = 0.0;
    protected int mMinIdentity = 50;

    protected int mPreFilterMinScore = 50;
    protected int mPreFilterMinOverlap = 50;
    protected int mPreFilterAlgorithm = -3;

    /**
     * @param filter the <code>OutputFilter</code>
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder outputFilter(final OutputFilter filter) {
      mOutputFilter = filter;
      return this;
    }

    /**
     * @param zip true iff output to be zipped.
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder zip(final boolean zip) {
      mZip = zip;
      return this;
    }

    /**
     * @param topN number of alternatives to be held.
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder topN(final int topN) {
      mTopN = topN;
      return this;
    }

    /**
     * @param maxTopResults number of results to output for topN.
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder maxTopResults(final int maxTopResults) {
      mMaxTopResults = maxTopResults;
      return this;
    }

    /**
     * @param exclude true iff to exclude repeated reads.
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder exclude(final boolean exclude) {
      mExcludeRepeats = exclude;
      return this;
    }

    /**
     * @param useIds true iff to use numeric identifiers for sequences instead of names.
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder useids(final boolean useIds) {
      mUseSequenceIds = useIds;
      return this;
    }

    /**
     * @param errorLimit maximum mismatch error allowed when outputting hits.
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder errorLimit(final int errorLimit) {
      if (errorLimit < 0 || errorLimit > MapFlags.MAX_SCORE) {
        throw new IllegalArgumentException(" errorLimit=" + errorLimit);
      }
      mErrorLimit = errorLimit;
      return this;
    }

    /**
     * @param maxMismatches unmated maximum mismatches allowed when outputting unmated hits.
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder unmatedMaxMismatches(final IntegerOrPercentage maxMismatches) {
      mUnmatedMaxMismatches = maxMismatches;
      return this;
    }

    /**
     * @param maxMismatches mated maximum mismatches allowed when outputting unmated hits.
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder matedMaxMismatches(final IntegerOrPercentage maxMismatches) {
      mMatedMaxMismatches = maxMismatches;
      return this;
    }

    /**
     * @param minIdentity minimum percent identity for output
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder minIdentity(final int minIdentity) {
      mMinIdentity = minIdentity;
      return this;
    }

    /**
     * @param minScore minimum percent identity in pre-filter
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder preFilterMinScore(final int minScore) {
      mPreFilterMinScore = minScore;
      return this;
    }

    /**
     * @param minOverlap minimum percent of the read that must overlap the template
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder preFilterMinOverlap(final int minOverlap) {
      mPreFilterMinOverlap = minOverlap;
      return this;
    }

    /**
     * @param algorithm algorithm used by pre-filter.
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder preFilterAlgorithm(final int algorithm) {
      mPreFilterAlgorithm = algorithm;
      return this;
    }

    /**
     * @param maxEScore maximum E-Score for output
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder maxEScore(final double maxEScore) {
      mMaxEScore = maxEScore;
      return this;
    }

    /**
     * @param minBitScore minimum bit score for output
     * @return this builder, so calls can be chained.
     */
    public NgsFilterParamsBuilder minBitScore(final double minBitScore) {
      mMinBitScore = minBitScore;
      return this;
    }

    /**
     * Creates a NgsParams using the current builder
     * configuration.
     * @return the new NgsParams
     */
    public NgsFilterParams create() {
      return new NgsFilterParams(this);
    }
  }



  private final OutputFilter mOutputFilter;

  private final int mTopN;

  private final int mMaxTopResults;

  private final boolean mExcludeRepeats;

  private final boolean mUseSequenceIds;

  private final boolean mZip;

  private final int mErrorLimit;
  private final IntegerOrPercentage mUnmatedMaxMismatches;
  private final IntegerOrPercentage mMatedMaxMismatches;
  private final int mMinIdentity;
  private final int mPreFilterMinScore;
  private final int mPreFilterMinOverlap;
  private final int mPreFilterAlgorithm;
  private final double mMaxEScore;
  private final double mMinBitScore;

  /**
   * @param builder the builder object.
   */
  public NgsFilterParams(final NgsFilterParamsBuilder builder) {
    mOutputFilter = builder.mOutputFilter;
    mZip = builder.mZip;
    mTopN = builder.mTopN;
    mMaxTopResults = builder.mMaxTopResults;
    mExcludeRepeats = builder.mExcludeRepeats;
    mUseSequenceIds = builder.mUseSequenceIds;
    mErrorLimit = builder.mErrorLimit;
    mUnmatedMaxMismatches = builder.mUnmatedMaxMismatches;
    mMatedMaxMismatches = builder.mMatedMaxMismatches;
    mMinIdentity = builder.mMinIdentity;
    mPreFilterMinScore = builder.mPreFilterMinScore;
    mPreFilterMinOverlap = builder.mPreFilterMinOverlap;
    mPreFilterAlgorithm = builder.mPreFilterAlgorithm;
    mMaxEScore = builder.mMaxEScore;
    mMinBitScore = builder.mMinBitScore;
  }


  /**
   * Check if output files to be zipped.
   * @return true if output to be zipped.
   */
  public boolean zip() {
    return mZip;
  }

  /**
   * Returns the current output filter to be used
   * @return output filter
   */
  public OutputFilter outputFilter() {
    return mOutputFilter;
  }

  /**
   * Get the error limit used with score-indel. Used in output to limit what is written out.
   *
   * @return the error limit (&ge;0).
   */
  public int errorLimit() {
    return mErrorLimit;
  }

  /**
   * Get the unmated max mismatches, used in output to limit unmated output
   * @return unmated error limit
   */
  public IntegerOrPercentage unmatedMaxMismatches() {
    return mUnmatedMaxMismatches;
  }

  /**
   * Get the mated max mismatches, used in output to limit unmated output
   * @return unmated error limit
   */
  public IntegerOrPercentage matedMaxMismatches() {
    return mMatedMaxMismatches;
  }

  /**
   * Get the minimum percent identity, used in output
   * @return minimum percent identity
   */
  public int minIdentity() {
    return mMinIdentity;
  }

  /**
   * Get the minimum score (roughly percent identity), used by the hit pre-filter before alignment.
   * @return minimum score
   */
  public int preFilterMinScore() {
    return mPreFilterMinScore;
  }

  /**
   * Get the minimum read overlap percentage, used by the hit pre-filter before alignment.
   * @return minimum score
   */
  public int preFilterMinOverlap() {
    return mPreFilterMinOverlap;
  }

  /**
   * Get the minimum score (roughly percent identity), used by the hit pre-filter before alignment.
   * @return minimum score
   */
  public int preFilterAlgorithm() {
    return mPreFilterAlgorithm;
  }

  /**
   * Get the maximum e-score, used in output
   * @return maximum e-score
   */
  public double maxEScore() {
    return mMaxEScore;
  }

  /**
   * Get the minimum bit score, used in output
   * @return minimum bit score
   */
  public double minBitScore() {
    return mMinBitScore;
  }

  /**
   * Get the number of repeated hits to be remembered for each read.
   * @return the number of repeated hits to be remembered for each read.
   */
  public int topN() {
    return mTopN;
  }

  /**
   * Get the number of results to be output for topN.
   * @return the number of results to be output for topN.
   */
  public int maxTopResults() {
    return mMaxTopResults;
  }

  /**
   * Check if reads with repeated hits should be excluded.
   * @return true if reads with repeated hits should be excluded.
   */
  public boolean exclude() {
    return mExcludeRepeats;
  }

  /**
   * Check if numeric identifiers should be used to name sequences (rather than the names of the sequences
   * as specified in the input files).
   * @return true iff numeric identifiers should be used.
   */
  public boolean useids() {
    return mUseSequenceIds;
  }

  @Override
  public int hashCode() {
    return Utils.hash(new Object[] {this.mOutputFilter, this.mZip, this.mTopN, this.mMaxTopResults, this.mExcludeRepeats, this.mUseSequenceIds, this.mErrorLimit, this.mUnmatedMaxMismatches, this.mMatedMaxMismatches});
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final NgsFilterParams that = (NgsFilterParams) obj;
    return Utils.equals(
        new Object[] {this.mOutputFilter, this.mZip, this.mTopN, this.mMaxTopResults, this.mExcludeRepeats, this.mUseSequenceIds, this.mErrorLimit, this.mUnmatedMaxMismatches, this.mMatedMaxMismatches},
        new Object[] {that.mOutputFilter, that.mZip, that.mTopN, this.mMaxTopResults, that.mExcludeRepeats, that.mUseSequenceIds, that.mErrorLimit, this.mUnmatedMaxMismatches, this.mMatedMaxMismatches}
    );
  }

  @Override
  public void toString(final StringBuilder sb) {
    toString("", sb);
  }

  /**
   * Output with prefix so that more readable when used with <code>Task</code> output.
   *
   * @param prefix to be added to each line.
   * @param sb where the result is to be placed.
   */
  public void toString(final String prefix, final StringBuilder sb) {
    sb.append(prefix).append("NgsFilterParams").append(" filter=").append(mOutputFilter).append(" topN=").append(mTopN).append(" maxTopResults=").append(mMaxTopResults).append(" error limit=").append(mErrorLimit).append(" mated max mismatches=").append(mMatedMaxMismatches).append(" unmated max mismatches=").append(mUnmatedMaxMismatches).append(" exclude=").append(mExcludeRepeats).append(" use-ids=").append(mUseSequenceIds).append(" zip=").append(mZip);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mTopN > 0);
    Exam.assertTrue(mMaxTopResults > 0);
    Exam.assertTrue(mErrorLimit >= 0);
    Exam.assertTrue(0 <= mMinIdentity && mMinIdentity <= 100);
    Exam.assertTrue(0 <= mPreFilterMinScore && mPreFilterMinScore <= 100);
    Exam.assertTrue(0 <= mPreFilterMinOverlap && mPreFilterMinOverlap <= 100);
    return true;
  }

}
