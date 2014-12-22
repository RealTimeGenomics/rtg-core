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
package com.rtg.sam;

import java.io.File;
import java.util.Collection;

/**
 * Encapsulates a single collection of mapping output files (typically when
 * used as input for another command).
 *
 */
public abstract class SingleMappedParams extends MappedParams {

  /**
   * Builder class for SingleMappedParams
   * @param <B> the builder type
   */
  public abstract static class SingleMappedParamsBuilder<B extends SingleMappedParamsBuilder<B>> extends MappedParamsBuilder<B> {

    protected Collection<File> mMapped;
    protected int mIoThreads = 1;
    protected int mExecThreads = 1;

    /**
     * Sets the mapping input files.
     * @param mapped the files containing mappings.
     * @return this builder, so calls can be chained.
     */
    public B mapped(final Collection<File> mapped) {
      mMapped = mapped;
      return self();
    }

    /**
     * The number of additional threads to use for Sam file reading and parsing.
     * @param ioThreads number of additional threads to run for reading Sam files.
     *          Default is 0.
     * @return this builder, so calls can be chained.
     */
    public B ioThreads(final int ioThreads) {
      mIoThreads = ioThreads;
      return self();
    }

    /**
     * The number of additional threads to use for Sam file reading and parsing.
     * @param execThreads number of additional threads to run for reading Sam files.
     *          Default is 0.
     * @return this builder, so calls can be chained.
     */
    public B execThreads(final int execThreads) {
      mExecThreads = execThreads;
      return self();
    }

  }

  private final Collection<File> mMapped;
  private final int mIoThreads;
  private final int mExecThreads;

  protected SingleMappedParams(final SingleMappedParamsBuilder<?> builder) {
    super(builder);
    mMapped = builder.mMapped;
    mIoThreads = builder.mIoThreads;
    mExecThreads = builder.mExecThreads;
  }

  /**
   * Get the mapped reads (in SAM/BAM format).
   * @return the mapped reads (in SAM/BAM format).
   */
  public Collection<File> mapped() {
    return mMapped;
  }

  /**
   * Get the number of extra threads to use in the Sam File Reading Pool
   * (Default 0)
   * @return the number of threads to use
   */
  public int ioThreads() {
    return mIoThreads;
  }

  /**
   * Get the number of threads to use for execution.
   * @return the number of threads to use for execution (default 1).
   */
  public int execThreads() {
    return mExecThreads;
  }
}
