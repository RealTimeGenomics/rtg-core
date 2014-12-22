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
package com.rtg.util.diagnostic;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory.
 *
 */
public class AllTests extends TestSuite {

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.util.diagnostic");

    suite.addTestSuite(CliDiagnosticListenerTest.class);
    suite.addTest(DiagnosticEventTest.suite());
    suite.addTest(DiagnosticTest.suite());
    suite.addTest(ErrorEventTest.suite());
    suite.addTest(ErrorTypeTest.suite());
    suite.addTestSuite(ListenerTypeTest.class);
    suite.addTestSuite(OneShotTimerTest.class);
    suite.addTest(ParallelProgressTest.suite());
    suite.addTest(SlimExceptionTest.suite());
    suite.addTestSuite(SpyCounterTest.class);
    suite.addTestSuite(SpyHistogramTest.class);
    suite.addTestSuite(SpyTimerTest.class);
    suite.addTestSuite(SpyTest.class);
    suite.addTestSuite(TalkbackTest.class);
    suite.addTestSuite(TimerTest.class);
    suite.addTest(WarningEventTest.suite());
    suite.addTest(WarningTypeTest.suite());
    suite.addTest(InformationEventTest.suite());
    suite.addTest(InformationTypeTest.suite());
    suite.addTestSuite(NoTalkbackSlimExceptionTest.class);
    suite.addTestSuite(WarningsTest.class);
    return suite;
  }
}

