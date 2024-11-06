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
package com.rtg.assembler;

import java.io.File;
import java.util.List;

import com.rtg.launcher.ModuleParams;
import com.rtg.util.intervals.LongRange;

/**
 * params type for DeBruijnParams
 */
public class DeBruijnParams extends ModuleParams {
  private final List<File> mInputFiles;
  private final int mKmerSize;
  private final int mMinHashFrequency;
  private final boolean mUseStringKmers;
  private final double mMergeRatio;
  private final File mDirectory;
  private final LongRange mRegion;

  protected DeBruijnParams(Builder builder) {
    super(builder);
    mInputFiles = builder.mInputFiles;
    mKmerSize = builder.mKmerSize;
    mMinHashFrequency = builder.mMinHashFrequency;
    mUseStringKmers = builder.mUseStringKmers;
    mMergeRatio = builder.mMergeRatio;
    mDirectory = builder.mDirectory;
    mRegion = builder.mRegion;
  }
  /**
   * @return the list of sequence files for graph building
   */
  public List<File> inputFiles() {
    return mInputFiles;
  }
  /**
   * @return the Kmer size you should use when constructing a graph
   */
  public int kmerSize() {
    return mKmerSize;
  }
  /**
   * @return minimum hash frequency or <code>-1</code> if you should use the automatically calculated hash frequency threshold
   */
  public int minHashFrequency() {
    return mMinHashFrequency;
  }
  /**
   * @return true if you should use string based Kmers
   */
  public boolean useStringKmers() {
    return mUseStringKmers;
  }
  /**
   * @return avoid merging bubbles were the ratio of kmers on the branches is below this
   */
  public double mergeRatio() {
    return mMergeRatio;
  }
  /**
   * @return restriction region for reads
   */
  public LongRange region() {
    return mRegion;
  }
  @Override
  public File directory() {
    return mDirectory;
  }
  @Override
  public File file(String name) {
    return new File(mDirectory, name);
  }
  @Override
  public String toString() {
    return "DeBruijnParams"
        + " inputFiles=" + mInputFiles
        + " kmerSize=" + mKmerSize
        + " minHashFrequency=" + mMinHashFrequency
        + " useStringKmers=" + mUseStringKmers
        + " mergeRatio=" + mMergeRatio
        + " directory=" + mDirectory
        + " region=" + mRegion
        + "";
  }
  /**
   * @return a new builder for this params class
   */
  public static Builder builder() {
    return new Builder();
  }
  /** Builder class for <code>DeBruijnParams</code> */
  public static class Builder extends ModuleParams.ModuleParamsBuilder<Builder> {
    private List<File> mInputFiles;
    private int mKmerSize;
    private int mMinHashFrequency;
    private boolean mUseStringKmers;
    private double mMergeRatio;
    private File mDirectory;
    private LongRange mRegion = LongRange.NONE;
    /**
     * @param inputFiles the list of sequence files for graph building
     * @return this so calls can be chained
     */
    public Builder inputFiles(List<File> inputFiles) {
      mInputFiles = inputFiles;
      return this;
    }
    /**
     * @param kmerSize the Kmer size you should use when constructing a graph
     * @return this so calls can be chained
     */
    public Builder kmerSize(int kmerSize) {
      mKmerSize = kmerSize;
      return this;
    }
    /**
     * @param minHashFrequency minimum hash frequency or <code>-1</code> if you should use the automatically calculated hash frequency threshold
     * @return this so calls can be chained
     */
    public Builder minHashFrequency(int minHashFrequency) {
      mMinHashFrequency = minHashFrequency;
      return this;
    }
    /**
     * @param useStringKmers true if you should use string based Kmers
     * @return this so calls can be chained
     */
    public Builder useStringKmers(boolean useStringKmers) {
      mUseStringKmers = useStringKmers;
      return this;
    }
    /**
     * @param mergeRatio avoid merging bubbles were the ratio of kmers on the branches is below this
     * @return this so calls can be chained
     */
    public Builder mergeRatio(double mergeRatio) {
      mMergeRatio = mergeRatio;
      return this;
    }
    /**
     * @param region reads restriction region
     * @return this so calls can be chained
     */
    public Builder region(LongRange region) {
      mRegion = region;
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
      @Override
      protected Builder self() {
        return this;
      }
      /**
       * @return an initialised <code>DeBruijnParams</code>
       */
      public DeBruijnParams create() {
        return new DeBruijnParams(this);
      }
  }
}
