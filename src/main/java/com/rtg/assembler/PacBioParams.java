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
