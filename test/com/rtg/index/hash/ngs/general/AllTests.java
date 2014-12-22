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
package com.rtg.index.hash.ngs.general;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory.
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.index.hash.ngs.general");
    suite.addTestSuite(CombinatorTest.class);
    suite.addTest(ExtractAbstractTest.suite());
    suite.addTestSuite(ExtractReadTest.class);
    suite.addTestSuite(ExtractTemplateTest.class);
    suite.addTestSuite(MaskAnalysisTest.class);
    suite.addTestSuite(MaskIndelCountTest.class);
    suite.addTestSuite(MaskTest.class);
    suite.addTestSuite(SingleMaskTest.class);
    suite.addTestSuite(SkelTest.class);
    suite.addTestSuite(SkeletonTest.class);
    suite.addTestSuite(SkeletonAltTest.class);
    suite.addTestSuite(UtilTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

