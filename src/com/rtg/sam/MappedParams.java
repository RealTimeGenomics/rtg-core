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

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.OutputModuleParams;
import com.rtg.util.test.params.ParamsNoField;

/**
 * MappedParams represents a Params for a high-level slim module with
 * a genome which the input SAM records were mapped against and
 * a SamFilterParams controlling its filtering of SAM input records.
 */
public abstract class MappedParams extends OutputModuleParams {

  /**
   * A builder class for <code>SamFilteredModuleParams</code>.
   * @param <B> builder type
   */
  public abstract static class MappedParamsBuilder<B extends MappedParamsBuilder<B>> extends OutputModuleParamsBuilder<B> {

    protected SamFilterParams mFilterParams = SamFilterParams.builder().create();
    protected ISequenceParams mGenome;

    boolean mIgnoreIncompatibleSamHeaders = false;

    /**
     * SAM filtering parameters.
     * @param val filtering parameters
     * @return this builder, so calls can be chained.
     */
    public B filterParams(final SamFilterParams val) {
      mFilterParams = val;
      return self();
    }

    /**
     * Sets the genome reader.
     * @param params the parameters used for output.
     * @return this builder, so calls can be chained.
     */
    public B genome(final ISequenceParams params) {
      mGenome = params;
      return self();
    }

    /**
     * true to ignore incompatible headers when merging SAM records
     * @param val the value
     * @return this builder
     */
    public B ignoreIncompatibleSamHeaders(boolean val) {
      mIgnoreIncompatibleSamHeaders = val;
      return self();
    }
  }

  private final SamFilterParams mFilterParams;
  private final ISequenceParams mGenome;
  private final boolean mIgnoreIncompatibleSamHeaders;

  /**
   * @param builder the builder object.
   */
  protected MappedParams(MappedParamsBuilder<?> builder) {
    super(builder);
    mFilterParams = builder.mFilterParams;
    mGenome = builder.mGenome;
    mIgnoreIncompatibleSamHeaders = builder.mIgnoreIncompatibleSamHeaders;
  }

  /**
   * Get the filter params object
   * @return the params
   */
  public SamFilterParams filterParams() {
    return mFilterParams;
  }

  /**
   *  Get parameters for the reference genome.
   * @return parameters for the reference genome.
   */
  public ISequenceParams genome() {
    return mGenome;
  }

  /**
   * @return true to ignore incompatible headers when merging SAM records
   */
  public boolean ignoreIncompatibleSamHeaders() {
    return mIgnoreIncompatibleSamHeaders;
  }

  @Override
  @ParamsNoField
  //params should be invariant.
  public void close() throws IOException {
    if (mGenome != null) {
      mGenome.close();
    }
  }

  @Override
  @ParamsNoField
  public boolean closed() {
    return mGenome == null || mGenome.closed();
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    sb.append(mFilterParams.toString()).append(LS);
    if (mGenome != null) {
      sb.append("    ").append(mGenome.toString());
      sb.append(LS);
    }
    sb.append(super.toString());
    return sb.toString();
  }
}
