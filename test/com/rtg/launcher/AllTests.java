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

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.launcher");
    suite.addTestSuite(DummyStatisticsTest.class);
    suite.addTestSuite(DefaultReaderParamsTest.class);
    suite.addTestSuite(DummyCliTest.class);
    suite.addTestSuite(GlobalFlagsTest.class);
    suite.addTestSuite(NoStatisticsTest.class);
    suite.addTestSuite(ParamsCliTest.class);
    suite.addTestSuite(SequenceParamsTest.class);
    suite.addTestSuite(OutputModuleParamsTest.class);
    suite.addTestSuite(OutputParamsTest.class);
    suite.addTestSuite(ModuleParamsTest.class);
    suite.addTestSuite(ParamsTaskTest.class);
    suite.addTestSuite(BuildParamsTest.class);
    suite.addTestSuite(ReaderParamsTest.class);
    suite.addTestSuite(BuildCommonTest.class);
    suite.addTestSuite(LoggedCliTest.class);
    suite.addTestSuite(HashingRegionTest.class);
    suite.addTestSuite(CommandLineFilesTest.class);
    suite.addTestSuite(CommonFlagsTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

