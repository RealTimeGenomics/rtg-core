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
package com.rtg.assembler;

import java.io.File;
import java.util.Collections;
import java.util.List;

import com.rtg.launcher.ModuleParams;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.intervals.LongRange;

/**
 * params type for AssembleParams
 */
public class AssembleParams extends ModuleParams {
  private final File mGraph;
  private final int mWordSize;
  private final int mStepSize;
  private final int mKmerSize;
  private final List<File> mReads;
  private final List<File> mReads454;
  private final List<File> mReadsMatePair;
  private final int mMinHashFrequency;
  private final double mMergeRatio;
  private final IntegerOrPercentage mMaxMismatches;
  private final int mMinInsertSize;
  private final int mMaxInsertSize;
  private final boolean mAlignments;
  private final int mMinPathReads;
  private final int mMinReadCount;
  private final int mConsensusThreshold;
  private final int mNumberThreads;
  private final File mDirectory;
  private final LongRange mRegion;

  protected AssembleParams(Builder builder) {
    super(builder);
    mGraph = builder.mGraph;
    mWordSize = builder.mWordSize;
    mStepSize = builder.mStepSize;
    mKmerSize = builder.mKmerSize;
    mReads = builder.mReads;
    mReads454 = builder.mReads454;
    mReadsMatePair = builder.mReadsMatePair;
    mMinHashFrequency = builder.mMinHashFrequency;
    mMergeRatio = builder.mMergeRatio;
    mMaxMismatches = builder.mMaxMismatches;
    mMinInsertSize = builder.mMinInsertSize;
    mMaxInsertSize = builder.mMaxInsertSize;
    mAlignments = builder.mAlignments;
    mMinPathReads = builder.mMinPathReads;
    mMinReadCount = builder.mMinReadCount;
    mConsensusThreshold = builder.mConsensusThreshold;
    mNumberThreads = builder.mNumberThreads;
    mDirectory = builder.mDirectory;
    mRegion = builder.mRegion;
  }
  /**
   * @return Start with this graph skipping kmer build
   */
  public File graph() {
    return mGraph;
  }
  /**
   * @return word size to use when mapping
   */
  public int wordSize() {
    return mWordSize;
  }
  /**
   * @return step size to use when mapping
   */
  public int stepSize() {
    return mStepSize;
  }
  /**
   * @return the Kmer size you should use when constructing a graph
   */
  public int kmerSize() {
    return mKmerSize;
  }
  /**
   * @return list of input reads
   */
  public List<File> reads() {
    return mReads;
  }
  /**
   * @return list of 454 input reads
   */
  public List<File> reads454() {
    return mReads454;
  }
  /**
   * @return list of 454 input reads
   */
  public List<File> readsMatePair() {
    return mReadsMatePair;
  }
  /**
   * @return minimum hash frequency or <code>-1</code> if you should use the automatically calculated hash frequency threshold
   */
  public int minHashFrequency() {
    return mMinHashFrequency;
  }
  /**
   * @return avoid merging bubbles where the ratio of kmers on the branches is below this
   */
  public double mergeRatio() {
    return mMergeRatio;
  }
  /**
   * @return number of mismatching bases allowed in alignment
   */
  public IntegerOrPercentage maxMismatches() {
    return mMaxMismatches;
  }
  /**
   * @return minimum insert size between fragments
   */
  public int minInsertSize() {
    return mMinInsertSize;
  }
  /**
   * @return maximum insert size between fragments
   */
  public int maxInsertSize() {
    return mMaxInsertSize;
  }
  /**
   * @return include output of alignments
   */
  public boolean alignments() {
    return mAlignments;
  }
  /**
   * @return minimum number of reads before a path is believable
   */
  public int minPathReads() {
    return mMinPathReads;
  }
  /**
   * @return minimum number of reads to preserve a link in the graph
   */
  public int minReadCount() {
    return mMinReadCount;
  }
  /**
   * @return number of reads necessary to have confidence in a path
   */
  public int consensusThreshold() {
    return mConsensusThreshold;
  }
  /**
   * @return number of threads to use in processing
   */
  public int numberThreads() {
    return mNumberThreads;
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
    return "AssembleParams"
        + " graph=" + mGraph
        + " wordSize=" + mWordSize
        + " stepSize=" + mStepSize
        + " kmerSize=" + mKmerSize
        + " reads=" + mReads
        + " reads454=" + mReads454
        + " readsMatePair=" + mReadsMatePair
        + " minHashFrequency=" + mMinHashFrequency
        + " mergeRatio=" + mMergeRatio
        + " maxMismatches=" + mMaxMismatches
        + " minInsertSize=" + mMinInsertSize
        + " maxInsertSize=" + mMaxInsertSize
        + " alignments=" + mAlignments
        + " minPathReads=" + mMinPathReads
        + " minReadCount=" + mMinReadCount
        + " consensusThreshold=" + mConsensusThreshold
        + " numberThreads=" + mNumberThreads
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
  /** Builder class for <code>AssembleParams</code> */
  public static class Builder extends ModuleParams.ModuleParamsBuilder<Builder> {
    private File mGraph;
    private int mWordSize;
    private int mStepSize;
    private int mKmerSize;
    private List<File> mReads = Collections.emptyList();
    private List<File> mReads454 = Collections.emptyList();
    private List<File> mReadsMatePair = Collections.emptyList();
    private int mMinHashFrequency;
    private double mMergeRatio;
    private IntegerOrPercentage mMaxMismatches = new IntegerOrPercentage(0);
    private int mMinInsertSize = Integer.MAX_VALUE;
    private int mMaxInsertSize = Integer.MIN_VALUE;
    private boolean mAlignments;
    private int mMinPathReads =  -1;
    private int mMinReadCount =  -1;
    private int mConsensusThreshold;
    private int mNumberThreads = 1;
    private File mDirectory;
    private LongRange mRegion = LongRange.NONE;
    /**
     * @param graph Start with this graph skipping kmer build
     * @return this so calls can be chained
     */
    public Builder graph(File graph) {
      mGraph = graph;
      return this;
    }
    /**
     * @param wordSize word size to use when mapping
     * @return this so calls can be chained
     */
    public Builder wordSize(int wordSize) {
      mWordSize = wordSize;
      return this;
    }
    /**
     * @param stepSize step size to use when mapping
     * @return this so calls can be chained
     */
    public Builder stepSize(int stepSize) {
      mStepSize = stepSize;
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
     * @param reads list of input reads
     * @return this so calls can be chained
     */
    public Builder reads(List<File> reads) {
      mReads = reads;
      return this;
    }
    /**
     * @param reads454 list of 454 input reads
     * @return this so calls can be chained
     */
    public Builder reads454(List<File> reads454) {
      mReads454 = reads454;
      return this;
    }
    /**
     * @param readsMatePair list of mate pair input reads
     * @return this so calls can be chained
     */
    public Builder readsMatePair(List<File> readsMatePair) {
      mReadsMatePair = readsMatePair;
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
     * @param mergeRatio avoid merging bubbles where the ratio of kmers on the branches is below this
     * @return this so calls can be chained
     */
    public Builder mergeRatio(double mergeRatio) {
      mMergeRatio = mergeRatio;
      return this;
    }
    /**
     * @param maxMismatches number of mismatching bases allowed in alignment
     * @return this so calls can be chained
     */
    public Builder maxMismatches(IntegerOrPercentage maxMismatches) {
      mMaxMismatches = maxMismatches;
      return this;
    }
    /**
     * @param minInsertSize minimum insert size between fragments
     * @return this so calls can be chained
     */
    public Builder minInsertSize(int minInsertSize) {
      mMinInsertSize = minInsertSize;
      return this;
    }
    /**
     * @param maxInsertSize maximum insert size between fragments
     * @return this so calls can be chained
     */
    public Builder maxInsertSize(int maxInsertSize) {
      mMaxInsertSize = maxInsertSize;
      return this;
    }
    /**
     * @param alignments include output of alignments
     * @return this so calls can be chained
     */
    public Builder alignments(boolean alignments) {
      mAlignments = alignments;
      return this;
    }
    /**
     * @param minPathReads minimum number of reads before a path is believable
     * @return this so calls can be chained
     */
    public Builder minPathReads(int minPathReads) {
      mMinPathReads = minPathReads;
      return this;
    }
    /**
     * @param minReadCount minimum number of reads to preserve a link in the graph
     * @return this so calls can be chained
     */
    public Builder minReadCount(int minReadCount) {
      mMinReadCount = minReadCount;
      return this;
    }
    /**
     * @param consensusThreshold number of reads necessary to have confidence in a path
     * @return this so calls can be chained
     */
    public Builder consensusThreshold(int consensusThreshold) {
      mConsensusThreshold = consensusThreshold;
      return this;
    }
    /**
     * @param numberThreads number of threads to use in processing
     * @return this so calls can be chained
     */
    public Builder numberThreads(int numberThreads) {
      mNumberThreads = numberThreads;
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
       * @return an initialised <code>AssembleParams</code>
       */
      public AssembleParams create() {
        return new AssembleParams(this);
      }
  }
}
