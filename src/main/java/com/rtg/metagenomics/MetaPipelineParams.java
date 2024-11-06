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

import com.rtg.launcher.OutputModuleParams;
import com.rtg.metagenomics.MetagenomicsWrapperCli.Platform;

/**
 * Parameters object for metagenomic pipelines.
 */
public class MetaPipelineParams extends OutputModuleParams {

  /**
   * Creates a <code>MetaPipelineParams</code> builder.
   * @return the builder.
   */
  public static MetaPipelineParamsBuilder builder() {
    return new MetaPipelineParamsBuilder();
  }

  /**
   * A builder class for <code>MetaPipelineParams</code>.
   */
  public static class MetaPipelineParamsBuilder extends OutputModuleParamsBuilder<MetaPipelineParamsBuilder> {

    protected File mProteinSdf = null;
    protected File mSpeciesSdf = null;
    protected File mFilterSdf = null;
    protected File mInputFile = null;
    protected File mInputLeft = null;
    protected File mInputRight = null;
    protected Platform mInputPlatform = Platform.ILLUMINA;

    /**
     * Set protein SDF.
     * @param proteinSdf The protein SDF to set.
     * @return this builder, so calls can be chained.
     */
    public MetaPipelineParamsBuilder proteinSdf(File proteinSdf) {
      mProteinSdf = proteinSdf;
      return self();
    }

    /**
     * Set species SDF.
     * @param speciesSdf The species SDF to set.
     * @return this builder, so calls can be chained.
     */
    public MetaPipelineParamsBuilder speciesSdf(File speciesSdf) {
      mSpeciesSdf = speciesSdf;
      return self();
    }

    /**
     * Set filter SDF.
     * @param filterSdf The filter SDF to set.
     * @return this builder, so calls can be chained.
     */
    public MetaPipelineParamsBuilder filterSdf(File filterSdf) {
      mFilterSdf = filterSdf;
      return self();
    }

    /**
     * Set input file.
     * @param inputFile The input file to set.
     * @return this builder, so calls can be chained.
     */
    public MetaPipelineParamsBuilder inputFile(File inputFile) {
      mInputFile = inputFile;
      return self();
    }

    /**
     * Set left input file.
     * @param inputLeft The left input file to set.
     * @return this builder, so calls can be chained.
     */
    public MetaPipelineParamsBuilder inputLeft(File inputLeft) {
      mInputLeft = inputLeft;
      return self();
    }

    /**
     * Set right input file.
     * @param inputRight The right input file to set.
     * @return this builder, so calls can be chained.
     */
    public MetaPipelineParamsBuilder inputRight(File inputRight) {
      mInputRight = inputRight;
      return self();
    }

    /**
     * Set machine type for platform of input file.
     * @param inputPlatform The machine type for platform of input file to set.
     * @return this builder, so calls can be chained.
     */
    public MetaPipelineParamsBuilder inputPlatform(Platform inputPlatform) {
      mInputPlatform = inputPlatform;
      return self();
    }

    @Override
    protected MetaPipelineParamsBuilder self() {
      return this;
    }

    /**
     * Creates a MetaPipelineParams using the current builder
     * configuration.
     * @return the new MetaPipelineParams
     */
    public MetaPipelineParams create() {
      return new MetaPipelineParams(this);
    }
  }

  private final File mProteinSdf;
  private final File mSpeciesSdf;
  private final File mFilterSdf;
  private final File mInputFile;
  private final File mInputLeft;
  private final File mInputRight;
  private final Platform mInputPlatform;

  /**
   * @param builder the builder object.
   */
  protected MetaPipelineParams(MetaPipelineParamsBuilder builder) {
    super(builder);
    mProteinSdf = builder.mProteinSdf;
    mSpeciesSdf = builder.mSpeciesSdf;
    mFilterSdf = builder.mFilterSdf;
    mInputFile = builder.mInputFile;
    mInputLeft = builder.mInputLeft;
    mInputRight = builder.mInputRight;
    mInputPlatform = builder.mInputPlatform;
  }

  /**
   * Get protein SDF.
   * @return Returns the protein SDF.
   */
  public File proteinSdf() {
    return mProteinSdf;
  }

  /**
   * Get species SDF.
   * @return Returns the species SDF.
   */
  public File speciesSdf() {
    return mSpeciesSdf;
  }

  /**
   * Get filter SDF.
   * @return Returns the filter SDF.
   */
  public File filterSdf() {
    return mFilterSdf;
  }

  /**
   * Get input file.
   * @return Returns the input file.
   */
  public File inputFile() {
    return mInputFile;
  }

  /**
   * Get left input file.
   * @return Returns the left input file.
   */
  public File inputLeft() {
    return mInputLeft;
  }

  /**
   * Get right input file.
   * @return Returns the right input file.
   */
  public File inputRight() {
    return mInputRight;
  }

  /**
   * Get machine type for platform of input file.
   * @return Returns machine type for platform of input file.
   */
  public Platform inputPlatform() {
    return mInputPlatform;
  }

}
