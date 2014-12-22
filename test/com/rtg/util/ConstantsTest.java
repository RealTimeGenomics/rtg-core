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
package com.rtg.util;


import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests the corresponding class.
 *
 */
public class ConstantsTest extends TestCase {

  /**
   */
  public ConstantsTest(final String name) {
    super(name);
  }
  public static Test suite() {
    return new TestSuite(ConstantsTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public void test() {
    assertEquals(1000, Constants.MINIMUM_FILE_CHUNK_SIZE);
    assertEquals(400, Constants.MAX_OPEN_FILES);
    assertEquals(1000000000, Constants.MAX_FILE_SIZE);
    assertTrue(Constants.MAX_FILE_SIZE >= Constants.MINIMUM_FILE_CHUNK_SIZE);
  }
}

