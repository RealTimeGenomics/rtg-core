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

package com.rtg.variant.sv.discord.pattern;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.util.List;

import com.rtg.launcher.ModuleParams;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.test.params.ParamsNoField;
import com.rtg.variant.sv.discord.DiscordantToolCli;

/**
 *         Date: 14/03/12
 *         Time: 3:34 PM
 */
public final class BreakpointPatternParams extends ModuleParams {
  /** Default value for same distance */
  public static final int DEFAULT_SAME_DISTANCE = 50;
  /** Default value for fragment length */
  public static final int DEFAULT_FRAGMENT_LENGTH = 500;
  private final  RegionRestriction mRegion;
  private final List<File> mFiles;
  private final File mDirectory;
  private final int mFragmentLength;
  private final int mSameDistance;

  private final int mMinDepth;

  private BreakpointPatternParams(Builder b) {
    super(b);
    mFiles = b.mFiles;
    mRegion = b.mRegion;
    mDirectory = b.mDirectory;
    mFragmentLength = b.mFragmentLength;
    mSameDistance = b.mSameDistance;
    mMinDepth = b.mMinDepth;
  }

  /**
   * @return a new builder object
   */
  public static Builder builder() {
    return new Builder();
  }

  /**
   * @return the region to operate on
   */
  public RegionRestriction region() {
    return mRegion;
  }

  /**
   * @return the VCF files to read
   */
  public List<File> files() {
    return mFiles;
  }

  /**
   * @return the maximum distance between breakpoints that will be considered the "same place"
   */
  public int sameDistance() {
    return mSameDistance;
  }
  /**
   * @return the maximum distance you'd expect to see the other end of a fragment
   */
  public int fragmentLength() {
    return mFragmentLength;
  }

  /**
   * @return the smallest number of reads contributing to a breakpoint for you to take it seriously
   */
  public int minDepth() {
    return mMinDepth;
  }


  @Override
  public File directory() {
    return mDirectory;
  }

  @Override
  @ParamsNoField
  public File file(String name) {
    return new File(mDirectory, name);
  }


  @Override
  @ParamsNoField
  public boolean closed() {
    return false;
  }

  @Override
  @ParamsNoField
  public void close() {
  }

  @Override
  public String toString() {
    return "BreakpointPatternParams"
      + " fragment-length=" + mFragmentLength
      + " same-distance=" + mSameDistance
      + " directory=" + mDirectory
      + " min-depth=" + mMinDepth
      + LS;
  }

  /**
   * Builder class for <code>BreakpointPatternParams</code>
   */
  public static final class Builder extends ModuleParamsBuilder<Builder> {

    private RegionRestriction mRegion;
    private List<File> mFiles;
    private File mDirectory;
    private int mSameDistance = DEFAULT_SAME_DISTANCE;
    private int mFragmentLength = DEFAULT_FRAGMENT_LENGTH;
    private int mMinDepth = DiscordantToolCli.DEFAULT_MIN_DEPTH;

    /**
     * @param files the VCF files to examine for patterns
     * @return this so calls can be chained
     */
    public Builder files(List<File> files) {
      mFiles = files;
      return this;
    }

    /**
     * @param region the region to examine for patterns
     * @return this so calls can be chained
     */
    public Builder region(RegionRestriction region) {
      mRegion = region;
      return this;
    }

    /**
     * @param distance the maximum distance between breakpoints that will be considered the "same place"
     * @return this so calls can be chained
     */
    public Builder sameDistance(int distance) {
      mSameDistance = distance;
      return this;
    }

    /**
     * @param length the maximum distance you'd expect to see the other end of a fragment
     * @return this so calls can be chained
     */
    public Builder fragmentLength(int length) {
      mFragmentLength = length;
      return this;
    }

    /**
     * @param directory output directory
     * @return this so calls can be chained
     */
    public Builder directory(File directory) {
      mDirectory = directory;
      return this;
    }

    /**
     * @param minDepth the smallest number of reads contributing to a breakpoint for you to take it seriously
     * @return this so calls can be chained
     */
    public Builder setMinDepth(int minDepth) {
      mMinDepth = minDepth;
      return this;
    }
    @Override
    protected Builder self() {
      return this;
    }

    /**
     * Creates a <code>BreakPointPatternParams</code> using the current builder
     * configuration.
     * @return the new <code>BreakPointPatternParams</code>
     */
    public BreakpointPatternParams create() {
      return new BreakpointPatternParams(this);
    }
  }
}
