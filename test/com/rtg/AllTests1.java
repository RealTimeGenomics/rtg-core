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

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * This class runs the first half of the SLIM tests (the a-o packages)
 * plus the ones in this package.
 *
 */
public class AllTests1 extends TestSuite {

  /**
   * @return <code>TestSuite</code>
   */
  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg");
    suite.addTestSuite(BlatantlyLeaveFileOpenTest.class);
    suite.addTestSuite(CommandLookupTest.class);
    suite.addTestSuite(SlimTest.class);
    suite.addTestSuite(ToolsEntryTest.class);
    suite.addTestSuite(CommandTest.class);
    suite.addTestSuite(CoreCommandTest.class);
    suite.addTestSuite(ToolsCommandTest.class);
    suite.addTestSuite(CommandCategoryTest.class);
    suite.addTestSuite(HelpCommandTest.class);
    suite.addTestSuite(ReleaseLevelTest.class);
    suite.addTestSuite(LicenseCommandTest.class);
    suite.addTestSuite(VersionCommandTest.class);
    // NOTE: keep these suites in alphabetical order.
    suite.addTest(com.rtg.alignment.AllTests.suite());
    suite.addTest(com.rtg.assembler.AllTests.suite());
    suite.addTest(com.rtg.bed.AllTests.suite());
    suite.addTest(com.rtg.calibrate.AllTests.suite());
    suite.addTest(com.rtg.complexity.AllTests.suite());
    suite.addTest(com.rtg.index.AllTests.suite());
    suite.addTest(com.rtg.launcher.AllTests.suite());
    suite.addTest(com.rtg.metagenomics.AllTests.suite());
    suite.addTest(com.rtg.misc.AllTests.suite());
    suite.addTest(com.rtg.ml.AllTests.suite());
    suite.addTest(com.rtg.mode.AllTests.suite());
    suite.addTest(com.rtg.ngs.AllTests.suite());
    suite.addTest(com.rtg.pairedend.AllTests.suite());
    suite.addTest(com.rtg.similarity.AllTests.suite());

    // Uncomment the following lines if you want to see when these tests finish.
    //    suite.addTest(new Test() {
    //
    //      @Override
    //      public int countTestCases() {
    //        return 0;
    //      }
    //
    //      @Override
    //      public void run(TestResult arg0) {
    //        System.out.println("Done1");
    //      }});
    return suite;
  }


  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
