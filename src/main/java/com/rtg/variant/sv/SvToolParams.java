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
package com.rtg.variant.sv;



import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.test.params.ParamsNoField;

/**
 */
public final class SvToolParams extends SvParams {

  static final String NAME_SIMPLE = "sv_simple.tsv";
  static final String NAME_BAYESIAN = "sv_bayesian.tsv";

  /**
   * Creates a <code>SvToolParamsBuilder</code>.
   * @return the builder.
   */
  public static SvToolParamsBuilder builder() {
    return new SvToolParamsBuilder();
  }

  /**
   * A builder class for <code>SvToolParams</code>.
   */
  public static final class SvToolParamsBuilder extends SvParamsBuilder<SvToolParamsBuilder> {

    boolean mOutputSimple;
    boolean mHeterozygous;

    int mBinSize;
    int mStepSize;
    int mFineStepSize;

    File mCorrectionsFile = null;

    @Override
    protected SvToolParamsBuilder self() {
      return this;
    }

    /**
     * Sets the bin size for printing simple signals (but not the main Bayesian signals).
     *
     * @param num number of nucleotides to put into one bin.
     * @return this builder, so calls can be chained.
     */
    public SvToolParamsBuilder binSize(final int num) {
      mBinSize = num;
      return self();
    }

    /**
     * Sets the step size
     * @param stepSize the step size
     * @return this builder, so calls can be chained.
     */
    public SvToolParamsBuilder stepSize(int stepSize) {
      mStepSize = stepSize;
      return self();
    }

    /**
     * Sets the step size to use in interesting regions
     * @param fineStepSize the step size for interesting regions
     * @return this builder, so calls can be chained.
     */
    public SvToolParamsBuilder fineStepSize(int fineStepSize) {
      mFineStepSize = fineStepSize;
      return self();
    }

    /**
     * Sets the corrections file.
     * @param file the with a correction value at each position.
     * @return this builder, so calls can be chained.
     */
    public SvToolParamsBuilder correctionsFile(File file) {
      mCorrectionsFile = file;
      return self();
    }

    /**
     * Set whether to output the simple signals.
     * @param simple true if simple output should be produced
     * @return this builder, so calls can be chained.
     */
    public SvToolParamsBuilder outputSimple(final boolean simple) {
      mOutputSimple = simple;
      return self();
    }

    /**
     * Set whether to also use heterozygous models.
     * @param hetero true if heterozygous models should be added
     * @return this builder, so calls can be chained.
     */
    public SvToolParamsBuilder heterozygous(final boolean hetero) {
      mHeterozygous = hetero;
      return self();
    }

    /**
     * Creates a <code>SVToolsParams</code> using the current builder
     * configuration.
     * @return the new <code>SVToolsParams</code>
     */
    public SvToolParams create() {
      return new SvToolParams(this);
    }

  }

  private final boolean mOutputSimple;
  private final boolean mHeterozygous;
  private final int mBinSize;
  private final int mStepSize;
  private final int mFineStepSize;

  private final File mCorrectionsFile;


  /**
   * @param builder the builder object.
   */
  SvToolParams(final SvToolParamsBuilder builder) {
    super(builder);
    mBinSize = builder.mBinSize;
    mStepSize = builder.mStepSize;
    mFineStepSize = builder.mFineStepSize;
    mOutputSimple = builder.mOutputSimple;
    mHeterozygous = builder.mHeterozygous;
    mCorrectionsFile = builder.mCorrectionsFile;
  }

  /**
   *  Get bin size.
   * @return bin size.
   */
  public int binSize() {
    return mBinSize;
  }

  /**
   * @return step size
   */
  public int stepSize() {
    return mStepSize;
  }
  /**
   * @return step size in interesting regions
   */
  public int fineStepSize() {
    return mFineStepSize;
  }

  /**
   * @return the stream for writing the bayesian signals.
   * @throws IOException whenever.
   */
  @ParamsNoField
  public OutputStream bayesianStream() throws IOException {
    return outStream(NAME_BAYESIAN);
  }

  /**
   * Get whether to produce simple signal output.
   * @return true if simple signal output should be produced.
   */
  public boolean outputSimple() {
    return mOutputSimple;
  }

  /**
   * Get whether heterozygous models should be used.
   * @return true if heterozygous models should be added
   */
  public boolean heterozygous() {
    return mHeterozygous;
  }

  /**
   * Get file with corrections.
   * @return file with corrections.
   */
  public File correctionsFile() {
    return mCorrectionsFile;
  }

  @Override
  public String toString() {
    final String pref = "    ";
    final StringBuilder sb = new StringBuilder();
    sb.append("SvToolParams" + " mapped-reads=").append(mapped()).append(" output-simple=").append(mOutputSimple).append(" heterozygous=").append(mHeterozygous).append(" bin-size=").append(mBinSize).append(" step-size=").append(mStepSize).append(" fine-step-size=").append(mFineStepSize).append(" correctionsFile=").append(mCorrectionsFile).append(LS);
    sb.append(filterParams()).append(LS);
    if (genome() != null) {
      sb.append(pref).append(genome());
      sb.append(LS);
    }
    sb.append(outputParams());
    sb.append(LS);
    return sb.toString();
  }
}
