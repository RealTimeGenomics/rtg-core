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

package com.rtg.launcher;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.test.params.ParamsNoField;

/**
 * OutputModuleParams represents a Params for a high-level slim module with an OutputParams controlling its output directory.
 */
public abstract class OutputModuleParams extends ModuleParams {

  /**
   * A builder class for <code>OutputModuleParams</code>.
   * @param <B> builder type
   */
  public abstract static class OutputModuleParamsBuilder<B extends OutputModuleParamsBuilder<B>> extends ModuleParamsBuilder<B> {

    protected OutputParams mOutputParams;

    /**
     * Sets the output parameters.
     *
     * @param params the parameters used for output.
     * @return this builder, so calls can be chained.
     */
    public B outputParams(final OutputParams params) {
      mOutputParams = params;
      return self();
    }

  }

  private final OutputParams mOutputParams;

  /**
   * @param builder the builder object.
   */
  protected OutputModuleParams(OutputModuleParamsBuilder<?> builder) {
    super(builder);
    mOutputParams = builder.mOutputParams;
  }

  /**
   * Get the output params object
   * @return the params
   */
  public OutputParams outputParams() {
    return mOutputParams;
  }

  @Override
  @ParamsNoField
  public File file(final String name) {
    return mOutputParams.file(name);
  }

  @Override
  @ParamsNoField
  public File directory() {
    return mOutputParams == null ? null : mOutputParams.directory();
  }

  /**
   * Return the file that will be used for output given a base name.
   * This obeys any settings for whether compression should be applied.
   * @param name file name
   * @return file to be used for output
   */
  @ParamsNoField
  public File outFile(String name) {
    return mOutputParams.outFile(name);
  }

  /**
   * Get a stream to an output file in the output directory.
   * This obeys any settings for whether compression should be applied.
   * @param name file name
   * @return the stream.
   * @throws IOException if an I/O error occurs.
   */
  @ParamsNoField
  public OutputStream outStream(String name) throws IOException {
    return mOutputParams.outStream(name);
  }

  /**
   * @return true if output is to be compressed in blocked compressed format
   */
  @ParamsNoField
  public boolean blockCompressed() {
    return mOutputParams.isBlockCompressed();
  }

  @Override
  public String toString() {
    return mOutputParams.toString() + LS;
  }

}
