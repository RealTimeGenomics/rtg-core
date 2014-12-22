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
package com.rtg.metagenomics;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.sam.SingleMappedParams;
import com.rtg.util.test.params.ParamsNoField;

/**
 */
public final class SpeciesParams extends SingleMappedParams {

  /**
   * @return a builder for species params.
   */
  public static SpeciesParamsBuilder builder() {
    return new SpeciesParamsBuilder();
  }

  /**
   * Builder for <code>SpeciesParams</code>
   */
  public static final class SpeciesParamsBuilder extends SingleMappedParamsBuilder<SpeciesParamsBuilder> {
    File mReferenceMap;
    int mMinIter;
    boolean mVerbose;
    boolean mPrintAll;
    double mMinConfidence;

    @Override
    protected SpeciesParamsBuilder self() {
      return this;
    }

    /**
     * Sets the minimum number of iterations to perform
     * @param minIter the minimum number of iterations to perform
     * @return this builder, so calls can be chained.
     */
    public SpeciesParamsBuilder minIter(int minIter) {
      mMinIter = minIter;
      return self();
    }

    /**
     * Turn verbose output on or off
     * @param verbose true if verbose output should be produced
     * @return this builder, so calls can be chained.
     */
    public SpeciesParamsBuilder verbose(boolean verbose) {
      mVerbose = verbose;
      return self();
    }

    /**
     * Set the file containing the reference to species name mapping file
     * @param referenceMap the file containing the reference to species name mapping file
     * @return this builder, so calls can be chained.
     */
    public SpeciesParamsBuilder referenceMap(File referenceMap) {
      mReferenceMap = referenceMap;
      return self();
    }

    /**
     * @param printAll true, if print non present species in output
     * @return this builder, so calls can be chained.
     */
    public SpeciesParamsBuilder printAll(boolean printAll) {
      mPrintAll = printAll;
      return self();
    }

    /**
     * @param minConfidence report species above this confidence score
     * @return this builder, so calls can be chained.
     */
    public SpeciesParamsBuilder minConfidence(double minConfidence) {
      mMinConfidence = minConfidence;
      return self();
    }

    /**
     * Creates a <code>SpeciesParams</code> using the current builder
     * configuration.
     * @return the new <code>SpeciesParams</code>
     */
    public SpeciesParams create() {
      return new SpeciesParams(this);
    }
  }

  private final int mMinIter;
  private final boolean mVerbose;
  private final File mReferenceMap;
  private final boolean mPrintAll;
  private final double mMinConfidence;

  private SpeciesParams(final SpeciesParamsBuilder builder) {
    super(builder);
    mMinIter = builder.mMinIter;
    mVerbose = builder.mVerbose;
    mReferenceMap = builder.mReferenceMap;
    mPrintAll = builder.mPrintAll;
    mMinConfidence = builder.mMinConfidence;
  }

  /**
   * @return the stream for writing the species TSV.
   * @throws IOException whenever.
   */
  @ParamsNoField
  public OutputStream speciesStream() throws IOException {
    return outStream("species.tsv");
  }

  /**
   * Get minimum number of iterations.
   * @return minimum number of iterations.
   */
  public int minIter() {
    return mMinIter;
  }

  /**
   * Get status of verbose flag.
   * @return true if verbose option to be on.
   */
  public boolean verbose() {
    return mVerbose;
  }

  /**
   * Get the file containing the reference to species name mapping file.
   * @return the file containing the reference to species name mapping.
   */
  public File referenceMap() {
    return mReferenceMap;
  }

  /**
   * @return print non present species as well;
   */
  public boolean printAll() {
    return mPrintAll;
  }

  /**
   * @return minimum confidence score for reporting
   */
  public double minConfidence() {
    return mMinConfidence;
  }

  @Override
  public String toString() {
    return "SpeciesParams " + super.toString();
  }
}
