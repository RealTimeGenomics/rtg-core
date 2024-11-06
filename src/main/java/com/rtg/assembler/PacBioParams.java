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

/**
 * params type for PacBioParams
 */
public class PacBioParams extends ModuleParams {
  private final List<File> mReads;
  private final File mGraph;
  private final boolean mTrimGraph;
  private final File mDirectory;

  protected PacBioParams(Builder builder) {
    super(builder);
    mReads = builder.mReads;
    mGraph = builder.mGraph;
    mTrimGraph = builder.mTrimGraph;
    mDirectory = builder.mDirectory;
  }
  /**
    * @return read input sdf
    */
  public List<File> reads() {
    return mReads;
  }
  /**
    * @return graph directory
    */
  public File graph() {
    return mGraph;
  }
  /**
    * @return should short disconnected contigs be trimmed
    */
  public boolean trimGraph() {
    return mTrimGraph;
  }
  /**
    * @return output directory
    */
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
    return "PacBioParams"
        + " reads=" + mReads
        + " graph=" + mGraph
        + " trimGraph=" + mTrimGraph
        + " directory=" + mDirectory
        + "";
  }
  /**
   * @return a new builder for this params class
   */
  public static Builder builder() {
    return new Builder();
  }
  /** Builder class for <code>PacBioParams</code> */
  public static class Builder extends ModuleParams.ModuleParamsBuilder<Builder> {
    private List<File> mReads;
    private File mGraph;
    private boolean mTrimGraph;
    private File mDirectory;
    /**
      * @param reads read input sdf
      * @return this so calls can be chained
      */
    public Builder reads(List<File> reads) {
      mReads = reads;
      return this;
    }
    /**
      * @param graph graph directory
      * @return this so calls can be chained
      */
    public Builder graph(File graph) {
      mGraph = graph;
      return this;
    }
    /**
      * @param trimGraph should short disconnected contigs be trimmed
      * @return this so calls can be chained
      */
    public Builder trimGraph(boolean trimGraph) {
      mTrimGraph = trimGraph;
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
       * @return an initialised <code>PacBioParams</code>
       */
      public PacBioParams create() {
        return new PacBioParams(this);
      }
  }
}
