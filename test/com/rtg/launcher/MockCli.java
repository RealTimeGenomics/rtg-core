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


import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.IORunnable;


/**
 */
public class MockCli extends ParamsCli<MockCliParams> {

  private final File mOutDir;

  /**
   * Create a mock command line with specified output directory.
   * @param outDir output directory
   */
  public MockCli(final File outDir) {
    mOutDir = outDir;
  }

  /**
   * Create a mock command line with no output directory.
   */
  public MockCli() {
    this(null);
  }

  @Override
  public String applicationName() {
    return "MockCli";
  }

  @Override
  protected void initFlags() {
    MockCliParams.makeFlags(mFlags);
  }

  @Override
  protected MockCliParams makeParams() {
    return new MockCliParams(mFlags);
  }

  @Override
  protected IORunnable task(MockCliParams params, final OutputStream out) {
    try {
      out.write("Mock task did something".getBytes());
    } catch (final IOException e) {
      // do nothing
    }
    return null;
  }

  @Override
  protected File outputDirectory() {
    if (mOutDir != null) {
      return mOutDir;
    }
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public String moduleName() {
    return "";
  }

}
