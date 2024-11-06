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
import java.util.Collections;
import java.util.List;

import com.rtg.assembler.graph.MutableGraph;
import com.rtg.launcher.ModuleParams;
import com.rtg.util.IntegerOrPercentage;

/**
 * params type for GraphMapParams
 */
public class GraphMapParams extends ModuleParams {
  private final int mWordSize;
  private final int mStepSize;
  private final MutableGraph mGraph;
  private final List<File> mReads;
  private final List<File> mReads454;
  private final List<File> mReadsMatePair;
  private final IntegerOrPercentage mMaxMismatches;
  private final int mMinInsertSize;
  private final int mMaxInsertSize;
  private final File mAlignmentFile;
  private final int mNumberThreads;
  private final File mDirectory;

  protected GraphMapParams(Builder builder) {
    super(builder);
    mWordSize = builder.mWordSize;
    mStepSize = builder.mStepSize;
    mGraph = builder.mGraph;
    mReads = builder.mReads;
    mReads454 = builder.mReads454;
    mReadsMatePair = builder.mReadsMatePair;
    mMaxMismatches = builder.mMaxMismatches;
    mMinInsertSize = builder.mMinInsertSize;
    mMaxInsertSize = builder.mMaxInsertSize;
    mAlignmentFile = builder.mAlignmentFile;
    mNumberThreads = builder.mNumberThreads;
    mDirectory = builder.mDirectory;
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
   * @return graph to map against
   */
  public MutableGraph graph() {
    return mGraph;
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
   * @return list of mate pair input reads
   */
  public List<File> readsMatePair() {
    return mReadsMatePair;
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
   * @return file to write alignments
   */
  public File alignmentFile() {
    return mAlignmentFile;
  }
  /**
   * @return number of threads to use in processing
   */
  public int numberThreads() {
    return mNumberThreads;
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
    return "GraphMapParams"
        + " wordSize=" + mWordSize
        + " stepSize=" + mStepSize
        + " graph=" + mGraph
        + " reads=" + mReads
        + " reads454=" + mReads454
        + " readsMatePair=" + mReadsMatePair
        + " maxMismatches=" + mMaxMismatches
        + " minInsertSize=" + mMinInsertSize
        + " maxInsertSize=" + mMaxInsertSize
        + " alignmentFile=" + mAlignmentFile
        + " numberThreads=" + mNumberThreads
        + " directory=" + mDirectory
        + "";
  }
  /**
   * @return a new builder for this params class
   */
  public static Builder builder() {
    return new Builder();
  }
  /** Builder class for <code>GraphMapParams</code> */
  public static class Builder extends ModuleParams.ModuleParamsBuilder<Builder> {
    private int mWordSize;
    private int mStepSize;
    private MutableGraph mGraph;
    private List<File> mReads = Collections.emptyList();
    private List<File> mReads454 = Collections.emptyList();
    private List<File> mReadsMatePair = Collections.emptyList();
    private IntegerOrPercentage mMaxMismatches = new IntegerOrPercentage(0);
    private int mMinInsertSize = Integer.MAX_VALUE;
    private int mMaxInsertSize = Integer.MIN_VALUE;
    private File mAlignmentFile;
    private int mNumberThreads = 1;
    private File mDirectory;
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
     * @param graph graph to map against
     * @return this so calls can be chained
     */
    public Builder graph(MutableGraph graph) {
      mGraph = graph;
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
     * @param alignmentFile file to write alignments
     * @return this so calls can be chained
     */
    public Builder alignmentFile(File alignmentFile) {
      mAlignmentFile = alignmentFile;
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
       * @return an initialised <code>GraphMapParams</code>
       */
      public GraphMapParams create() {
        return new GraphMapParams(this);
      }
  }
}
