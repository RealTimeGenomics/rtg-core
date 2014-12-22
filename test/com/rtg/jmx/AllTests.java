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
package com.rtg.jmx;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class AllTests extends TestCase {

  public static Test suite() {
    final TestSuite suite = new TestSuite();
    suite.addTestSuite(ExternalCommandTest.class);
    suite.addTestSuite(MonUtilsTest.class);

    suite.addTestSuite(DiskStatsTest.class);
    suite.addTestSuite(NetworkStatsTest.class);
    suite.addTestSuite(MBeanStatsTest.class);
    suite.addTestSuite(ProgressStatsTest.class);

    suite.addTestSuite(LocalStatsTest.class);
    suite.addTestSuite(RecordStatsTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
