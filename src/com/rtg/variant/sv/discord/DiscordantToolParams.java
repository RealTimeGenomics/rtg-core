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
    boolean mIntersectionOnly;
//    Double mMaxAmbiguity;
//    Integer mMaxCoverage;

    int mMinBreakpointDepth;

    @Override
    protected DiscordantToolParamsBuilder self() {
      return this;
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

//    /**
//     * Sets the maximum proportion of ambiguous reads allowed
//     *
//     * @param ambiguity highest ambiguity allowed in output
//     * @return this builder, so calls can be chained.
//     */
//    public DiscordantToolParamsBuilder maxAmbiguity(Double ambiguity) {
//      mMaxAmbiguity = ambiguity;
//      return self();
//    }
//    /**
//     * Sets the maximum coverage allowed at breakpoints
//     *
//     * @param coverage highest coverage allowed in output
//     * @return this builder, so calls can be chained.
//     */
//    public DiscordantToolParamsBuilder maxCoverage(Integer coverage) {
//      mMaxCoverage = coverage;
//      return self();
//    }

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
  private final boolean mIntersectionOnly;
//  private final Double mMaxAmbiguity;
//  private final Integer mMaxCoverage;

  private final int mMinBreakpointDepth;

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
//    mMaxAmbiguity = builder.mMaxAmbiguity;
//    mMaxCoverage = builder.mMaxCoverage;
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

  @Override
  public String toString() {
    final String pref = "    ";
    final StringBuilder sb = new StringBuilder();
    sb.append("DiscordantToolParams" + " mapped-reads=").append(mapped().toString()).append(" output-tabix-index=").append(mOutputTabixIndex).append(" debug-output=").append(mDebugOutput).append(" bed-output=").append(mBedOutput).append(" intersection-only=").append(mIntersectionOnly).append(" min-breakpoint-depth=").append(mMinBreakpointDepth).append(LS);
    sb.append(filterParams().toString()).append(LS);
    if (genome() != null) {
      sb.append(pref).append(genome().toString());
      sb.append(LS);
    }
    sb.append(outputParams().toString());
    sb.append(LS);
    return sb.toString();
  }
}
