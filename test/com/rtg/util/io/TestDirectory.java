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
package com.rtg.util.io;

import java.io.File;
import java.io.IOException;

import com.rtg.util.test.FileHelper;

import junit.framework.Assert;

/**
 * Class to make it easier to create temp folder to unit tests.
 * Automatically adds an appropriate name, and the close method deletes the folder
 * and all its contents
 */
public class TestDirectory extends File implements AutoCloseable {

  /**
   * Constructs with suffix using name from invoking method
   * @throws IOException it happens maybe
   */
  public TestDirectory() throws IOException {
    this(getTestName());
  }

  /**
   * @param suffix suffix for temp dir file name
   * @throws IOException it happens maybe
   */
  public TestDirectory(String suffix) throws IOException {
    this(FileUtils.createTempDir("unitTest", suffix));
  }

  /**
   * @param dir directory to use for tests
   */
  public TestDirectory(File dir) {
    super(dir.getPath());
  }

  /**
   * deletes temporary directory and all contents
   */
  @Override
  public void close() {
    if (System.getProperty("testdirectory.nocleanup") != null) {
      System.err.println("Not deleting test directory " + this.toString());
    } else {
      Assert.assertTrue(FileHelper.deleteAll(this));
    }
  }

  /**
   * for use by constructor only
   * @return class and method name of constructor invoker
   */
  private static String getTestName() {
    final StackTraceElement[] elements = (new Exception()).getStackTrace();
    if (elements.length > 2) {
      return elements[2].getClassName() + "." + elements[2].getMethodName();
    } else {
      return "Unknown";
    }
  }
}
