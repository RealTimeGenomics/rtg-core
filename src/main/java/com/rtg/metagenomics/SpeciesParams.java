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
    IdentifierCreator mIdentifierCreator;

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
     * @param identifierCreator class to create identifiers
     * @return this builder, so calls can be chained.
     */
    public SpeciesParamsBuilder identifierCreator(final IdentifierCreator identifierCreator) {
      mIdentifierCreator = identifierCreator;
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
  private final IdentifierCreator mIdentifierCreator;

  private SpeciesParams(final SpeciesParamsBuilder builder) {
    super(builder);
    mMinIter = builder.mMinIter;
    mVerbose = builder.mVerbose;
    mReferenceMap = builder.mReferenceMap;
    mPrintAll = builder.mPrintAll;
    mMinConfidence = builder.mMinConfidence;
    mIdentifierCreator = builder.mIdentifierCreator;
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

  /**
   * @return class for creating identifiers.
   */
  public IdentifierCreator identifierCreator() {
    return mIdentifierCreator;
  }

  @Override
  public String toString() {
    return "SpeciesParams " + super.toString();
  }
}
