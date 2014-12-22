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
package com.rtg;

import com.rtg.util.TestUtils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests corresponding class
 */
public class CommandCategoryTest extends TestCase {

  public void testClass() {
    TestUtils.testEnum(CommandCategory.class, "[FORMAT, MAPPING, PROTEIN, VARIANT, METAGENOMICS, PIPELINES, SIMULATE, UTILITY, ASSEMBLY]");
  }

  public static Test suite() {
    return new TestSuite(CommandCategoryTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(new TestSuite(CommandCategoryTest.class));
  }
}
