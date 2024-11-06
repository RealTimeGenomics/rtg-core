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
package com.rtg.variant.cnv;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.util.Collection;

import com.rtg.sam.MappedParams;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 * Params for CNV module which takes two sets of input SAM mapping records.
 */
public final class CnvProductParams extends MappedParams implements Integrity {

  /**
   * Creates a <code>CnvProductParams</code> builder.
   * @return the builder.
   */
  public static CnvProductParamsBuilder builder() {
    return new CnvProductParamsBuilder();
  }

  /**
   * A builder class for <code>CnvProductParams</code>.
   */
  public static final class CnvProductParamsBuilder extends MappedParamsBuilder<CnvProductParamsBuilder> {

    Collection<File> mMappedBase;

    Collection<File> mMappedTarget;

    int mBucketSize = 100;

    boolean mFilterStartPositions = true;

    double mMagicConstant = 3.0;

    double mMultiplicationFactor = 3.0;

    boolean mExtraPenaltyOff = false;

    double mDivisionFactor = 3.0;

    int mThreads = 1;

    @Override
    protected CnvProductParamsBuilder self() {
      return this;
    }

    /**
     * Sets the mapping input files for the baseline of the comparison
     *
     * @param mapped the file containing mappings.
     * @return this builder, so calls can be chained.
     */
    public CnvProductParamsBuilder mappedBase(final Collection<File> mapped) {
      mMappedBase = mapped;
      return self();
    }

    /**
     * Sets the mapping input files for the target of the comparison
     *
     * @param mapped the file containing mappings.
     * @return this builder, so calls can be chained.
     */
    public CnvProductParamsBuilder mappedTarget(final Collection<File> mapped) {
      mMappedTarget = mapped;
      return self();
    }

    /**
     * Sets the bucket size.
     *
     * @param size the bucket size.
     * @return this builder, so calls can be chained.
     */
    public CnvProductParamsBuilder bucketSize(final int size) {
      mBucketSize = size;
      return self();
    }

    /**
     * Whether alignments at same start position should be treated as duplicates.
     * @param filterStart as above
     * @return this builder, so calls can be chained.
     */
    public CnvProductParamsBuilder filterStartPositions(final boolean filterStart) {
      mFilterStartPositions = filterStart;
      return self();
    }

    /**
     * It's magic
     * @param val new magic value
     * @return this builder, so calls can be chained.
     */
    public CnvProductParamsBuilder magicConstant(final double val) {
      mMagicConstant = val;
      return self();
    }

    /**
     * The factor to multiply the mean by for filtering out values above
     * @param val new multiplication factor
     * @return this builder, so calls can be chained.
     */
    public CnvProductParamsBuilder multiplicationFactor(final double val) {
      mMultiplicationFactor = val;
      return self();
    }

    /**
     * Switch extra penalty off
     * @param flag switch it off if true
     * @return this builder, so calls can be chained.
     */
    public CnvProductParamsBuilder extraPenaltyOff(final boolean flag) {
      mExtraPenaltyOff = flag;
      return self();
    }

    /**
     * The factor to divide the mean by for filtering out values below
     * @param val new division factor
     * @return this builder, so calls can be chained.
     */
    public CnvProductParamsBuilder divisionFactor(final double val) {
      mDivisionFactor = val;
      return self();
    }

    /**
     * Sets the number of additional threads to use for Sam file reading / parsing.
     *
     * @param threads the number of threads to employ in the Sam file reading. The default is 0.
     * @return this builder, so calls can be chained.
     */
    public CnvProductParamsBuilder threads(final int threads) {
      mThreads = threads;
      return self();
    }

    /**
     * Creates a <code>CnvProductParams</code> using the current builder
     * configuration.
     * @return the new <code>CnvProductParams</code>
     */
    public CnvProductParams create() {
      return new CnvProductParams(this);
    }
  }


  private final Collection<File> mMappedBase;
  private final Collection<File> mMappedTarget;
  private final int mBucketSize;
  private final boolean mFilterStartPositions;
  private final double mMagicConstant;
  private final double mMultiplicationFactor;
  private final boolean mExtraPenaltyOff;
  private final double mDivisionFactor;
  private final int mThreads;

  /**
   * @param builder the builder object.
   */
  CnvProductParams(final CnvProductParamsBuilder builder) {
    super(builder);
    mBucketSize = builder.mBucketSize;
    mMappedBase = builder.mMappedBase;
    mMappedTarget = builder.mMappedTarget;
    mFilterStartPositions = builder.mFilterStartPositions;
    mMagicConstant = builder.mMagicConstant;
    mMultiplicationFactor = builder.mMultiplicationFactor;
    mDivisionFactor = builder.mDivisionFactor;
    mThreads = builder.mThreads;
    mExtraPenaltyOff = builder.mExtraPenaltyOff;
  }

  /**
   * Get parameters for the bucket size.
   * @return parameters for the bucket size.
   */
  public int bucketSize() {
    return mBucketSize;
  }

  /**
   * SAM files for alignments representing the baseline
   * @return the files
   */
  public Collection<File> mappedBase() {
    return mMappedBase;
  }

  /**
   * SAM files for alignments representing the target of the comparison.
   * @return the files
   */
  public Collection<File> mappedTarget() {
    return mMappedTarget;
  }

  /**
   * Whether alignments at same start position should be treated as duplicates.
   * @return true for yes, false for no
   */
  public boolean filterStartPositions() {
    return mFilterStartPositions;
  }

  /**
   * It's magic
   * @return magic constant
   */
  public double magicConstant() {
    return mMagicConstant;
  }

  /**
   * The factor to multiply the mean by for filtering out values above
   * @return multiplication factor
   */
  public double multiplicationFactor() {
    return mMultiplicationFactor;
  }

  /**
   * The flag to switch extra penalty off
   * @return flag if extra penalty is off
   */
  public boolean extraPenaltyOff() {
    return mExtraPenaltyOff;
  }

  /**
   * The factor to divide the mean by for filtering out values below
   * @return division factor
   */
  public double divisionFactor() {
    return mDivisionFactor;
  }

  /**
   * Get the number of threads to be used for Sam file reading / parsing.
   * @return the number of threads to be used.
   */
  public int threads() {
    return mThreads;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mappedBase());
    Exam.assertNotNull(mappedTarget());
    Exam.assertTrue(mThreads >= 0);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    Exam.assertTrue(integrity());
    return true;
  }

  @Override
  public String toString() {
    return "CnvProductParams"
      + " baseLineInput=" + mappedBase()
      + " targetInput=" + mappedTarget()
      + " bucketSize=" + mBucketSize
      + " filterStartPositions=" + mFilterStartPositions
      + " divisionFactor=" + mDivisionFactor
      + " multiplicationFactor=" + mMultiplicationFactor
      + LS + super.toString();
  }

}
