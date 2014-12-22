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
    return ("CnvProductParams"
      + " baseLineInput=" + mappedBase().toString()
      + " targetInput=" + mappedTarget().toString()
      + " bucketSize=" + mBucketSize
      + " filterStartPositions=" + mFilterStartPositions
      + " divisionFactor=" + mDivisionFactor
      + " multiplicationFactor=" + mMultiplicationFactor
      + LS) + super.toString();
  }

}
