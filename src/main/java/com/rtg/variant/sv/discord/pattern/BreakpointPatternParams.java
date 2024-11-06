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

package com.rtg.variant.sv.discord.pattern;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.util.List;

import com.rtg.launcher.ModuleParams;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.test.params.ParamsNoField;

/**
 */
public final class BreakpointPatternParams extends ModuleParams {
  /** Default value for same distance */
  static final int DEFAULT_SAME_DISTANCE = 50;
  /** Default value for fragment length */
  static final int DEFAULT_FRAGMENT_LENGTH = 500;
  /** Default value for minimum support */
  static final Integer DEFAULT_MIN_DEPTH = 3;

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
    private int mMinDepth = DEFAULT_MIN_DEPTH;

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
