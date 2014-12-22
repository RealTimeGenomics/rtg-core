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
package com.rtg.scheduler;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.scheduler");

    suite.addTest(com.rtg.scheduler.enumtime.AllTests.suite());
    suite.addTest(com.rtg.scheduler.example.AllTests.suite());

    suite.addTestSuite(ExecutorRandomTest.class);
    suite.addTestSuite(ExecutorSequentialTest.class);
    suite.addTestSuite(ExecutorThreadedTest.class);
    suite.addTestSuite(JobTest.class);
    suite.addTestSuite(LookAheadTest.class);
    suite.addTestSuite(ResultTest.class);
    suite.addTestSuite(SchedulerSynchronizedTest.class);
    suite.addTestSuite(UtilTest.class);
    return suite;
  }
}
