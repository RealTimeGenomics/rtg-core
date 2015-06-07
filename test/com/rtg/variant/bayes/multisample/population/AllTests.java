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

package com.rtg.variant.bayes.multisample.population;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite(AllTests.class.getName());
    suite.addTestSuite(AlleleCountsFileReaderTest.class);
    suite.addTestSuite(AlleleCountsTest.class);
    suite.addTestSuite(ConvergencePlotTest.class);
    suite.addTestSuite(ConvergenceTest.class);
    suite.addTestSuite(DescriptionCountsTest.class);
    suite.addTestSuite(EmAlgorithmTest.class);
    suite.addTestSuite(EmResultTest.class);
    suite.addTestSuite(HwEstimatorTest.class);
    suite.addTestSuite(NullEstimatorTest.class);
    suite.addTestSuite(PopulationCliTest.class);
    suite.addTestSuite(PopulationCallerTest.class);
    suite.addTestSuite(PopulationNanoTest.class);
    suite.addTestSuite(PopulationHwHypothesesCreatorTest.class);
    suite.addTestSuite(PopulationCallerConfigurationTest.class);
    return suite;
  }
}
