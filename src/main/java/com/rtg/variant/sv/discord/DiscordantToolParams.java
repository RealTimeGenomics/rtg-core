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

package com.rtg.variant.sv.discord;

import static com.rtg.util.StringUtils.LS;

import com.rtg.variant.sv.SvParams;

/**
 */
public final class DiscordantToolParams extends SvParams {

  /**
   * Creates a <code>SvToolParamsBuilder</code>.
   * @return the builder.
   */
  public static DiscordantToolParamsBuilder builder() {
    return new DiscordantToolParamsBuilder();
  }

  /**
   * A builder class for <code>SvToolParams</code>.
   */
  public static final class DiscordantToolParamsBuilder extends SvParamsBuilder<DiscordantToolParamsBuilder> {

    boolean mOutputTabixIndex = true;
    boolean mDebugOutput = false;
    boolean mBedOutput = false;
    boolean mMultisample = false;
    boolean mIntersectionOnly = false;
    DiscordantTool.BamType mBamOutput = DiscordantTool.BamType.NONE;
    int mMinBreakpointDepth;
    double mOverlapFraction = 0;
    double mNumDeviations = 4.0;

    @Override
    protected DiscordantToolParamsBuilder self() {
      return this;
    }

    /**
     * @param num the number of standard deviations delineating concordant vs discordant fragment lengths
     * @return this builder, so calls can be chained.
     */
    public DiscordantToolParamsBuilder numDeviations(double num) {
      mNumDeviations = num;
      return self();
    }

    /**
     * @param val true if TABIX index should be created for output files.
     * @return this builder, so calls can be chained.
     */
    public DiscordantToolParamsBuilder outputTabixIndex(boolean val) {
      mOutputTabixIndex = val;
      return self();
    }

    /**
     * Set whether default output is debug
     * @param debug true = debug output, false = VCF output
     * @return this builder, so calls can be chained.
     */
    public DiscordantToolParamsBuilder debugOutput(final boolean debug) {
      mDebugOutput = debug;
      return self();
    }

    /**
     * Set whether bed output is produced
     * @param bed if true BED will be output
     * @return this builder, so calls can be chained.
     */
    public DiscordantToolParamsBuilder bedOutput(final boolean bed) {
      mBedOutput = bed;
      return self();
    }

    /**
     * Assume this fraction of the read length may have been erroneously aligned across a breakpoint due to alignment penalties.
     * Actual number varies according to alignment penalties and propensity for soft-clipping the dirty ends of alignments.
     * @param overlapFraction the fraction of read length to assume may overlap a breakend
     * @return this builder, so calls can be chained.
     */
    public DiscordantToolParamsBuilder overlapFraction(double overlapFraction) {
      mOverlapFraction = overlapFraction;
      return self();
    }

    /**
     * Set the minimum read depth threshold for breakpoints
     * @param depth breakpoints with less than this number of reads will be suppressed
     * @return this builder, so calls can be chained.
     */
    public DiscordantToolParamsBuilder minBreakpointDepth(final int depth) {
      mMinBreakpointDepth = depth;
      return self();
    }

    /**
     * @param intersectionOnly true if only breakpoints with a non-zero intersection should be output
     * @return this builder, so calls can be chained.
     */
    public DiscordantToolParamsBuilder intersectionOnly(final boolean intersectionOnly) {
      mIntersectionOnly = intersectionOnly;
      return self();
    }

    /**
     * @param multisample true if data from multiple samples can be provided
     * @return this builder, so calls can be chained.
     */
    public DiscordantToolParamsBuilder allowMultisample(boolean multisample) {
      mMultisample = multisample;
      return this;
    }

    /**
     * @param output the type of BAM output to produce
     * @return this builder, so calls can be chained.
     */
    public DiscordantToolParamsBuilder bamOutput(DiscordantTool.BamType output) {
      mBamOutput = output;
      return this;
    }

    /**
     * Creates a <code>SVToolsParams</code> using the current builder
     * configuration.
     * @return the new <code>SVToolsParams</code>
     */
    public DiscordantToolParams create() {
      return new DiscordantToolParams(this);
    }
  }


  private final boolean mOutputTabixIndex;
  private final boolean mDebugOutput;
  private final boolean mBedOutput;
  private final boolean mMultisample;
  private final boolean mIntersectionOnly;
  private final DiscordantTool.BamType mBamOutput;
  private final int mMinBreakpointDepth;
  private final double mOverlapFraction;
  private final double mNumDeviations;


  /**
   * @param builder the builder object.
   */
  DiscordantToolParams(DiscordantToolParamsBuilder builder) {
    super(builder);
    mOutputTabixIndex = builder.mOutputTabixIndex;
    mDebugOutput = builder.mDebugOutput;
    mBedOutput = builder.mBedOutput;
    mIntersectionOnly = builder.mIntersectionOnly;
    mMinBreakpointDepth = builder.mMinBreakpointDepth;
    mOverlapFraction = builder.mOverlapFraction;
    mMultisample = builder.mMultisample;
    mBamOutput = builder.mBamOutput;
    mNumDeviations = builder.mNumDeviations;
  }

  /**
   * @return the number of standard deviations delineating concordant vs discordant fragment lengths
   */
  public double numDeviations() {
    return mNumDeviations;
  }

  /**
   * @return true if TABIX index should be output.
   */
  public boolean outputTabixIndex() {
    return mOutputTabixIndex;
  }

  /**
   * Get whether we going to produce debug output
   * @return true if a debug file will be produced
   */
  public boolean debugOutput() {
    return mDebugOutput;
  }

  /**
   * @return true if multiple samples are permitted
   */
  public boolean multisample() {
    return mMultisample;
  }

  /**
   * @return what type of BAM to output
   */
  public DiscordantTool.BamType bamOutput() {
    return mBamOutput;
  }

  /**
   * Get whether we going to produce BED output
   * @return true if a BED file will be produced
   */
  public boolean bedOutput() {
    return mBedOutput;
  }

  /**
   * @return true if only breakpoints with a non-zero intersection should be output
   */
  public boolean intersectionOnly() {
    return mIntersectionOnly;
  }

  /**
   * @return the minimum number of reads a breakpoint must have in order to be output
   */
  public int minBreakpointDepth() {
    return mMinBreakpointDepth;
  }

  /**
   * @return the fraction of the read length that may have been erroneously aligned across a breakpoint due to alignment penalties.
   */
  public double overlapFraction() {
    return mOverlapFraction;
  }

  @Override
  public String toString() {
    final String pref = "    ";
    final StringBuilder sb = new StringBuilder();
    sb.append("DiscordantToolParams")
      .append(" mapped-reads=").append(mapped())
      .append(" output-tabix-index=").append(mOutputTabixIndex)
      .append(" debug-output=").append(mDebugOutput)
      .append(" bed-output=").append(mBedOutput)
      .append(" bam-output=").append(mBamOutput)
      .append(" intersection-only=").append(mIntersectionOnly)
      .append(" min-breakpoint-depth=").append(mMinBreakpointDepth)
      .append(" overlap-fraction=").append(mOverlapFraction).append(LS);
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
