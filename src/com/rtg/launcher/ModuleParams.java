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

import com.rtg.util.ObjectParams;
import com.rtg.util.ParamsBuilder;
import com.rtg.util.test.params.ParamsNoField;

/**
 * ModuleParams represents a Params for a high-level slim module. Currently these modules also need to have some notion of a main output directory (this should be abstracted out).
 */
public abstract class ModuleParams extends ObjectParams implements OutputDirParams {

  /**
   * @param <B> the builder type
   */
  public abstract static class ModuleParamsBuilder<B extends ModuleParamsBuilder<B>> extends ParamsBuilder<B> {

    protected String mName = "ModuleParams";

    /**
     * Sets the application name.
     * @param name the application name.
     * @return this builder, so calls can be chained.
     */
    public B name(final String name) {
      mName = name;
      return self();
    }

  }

  private final String mName;

  protected ModuleParams(final String name) {
    if (name == null) {
      throw new NullPointerException();
    }
    mName = name;
    mObjects = new Object[] {mName};
  }

  protected ModuleParams(final ModuleParamsBuilder<?> builder) {
    mName = builder.mName;
    mObjects = new Object[] {mName};
  }

  /**
   * Generate a child file in the output directory.
   * @param name of the child.
   * @return a child file in the output directory.
   */
  @ParamsNoField
  public abstract File file(final String name);

  /**
   * @return the external name of the application.
   */
  public String name() {
    return mName;
  }

  @Override
  @ParamsNoField
  public void close() throws IOException {
    // do nothing
  }

  @Override
  public String toString() {
    return name() + LS;
  }

}
