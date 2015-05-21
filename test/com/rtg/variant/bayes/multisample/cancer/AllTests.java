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
package com.rtg.variant.bayes.multisample.cancer;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite(AllTests.class.getName());

    suite.addTestSuite(ContaminatedSomaticCallerTest.class);
    suite.addTestSuite(PureSomaticCallerTest.class);
    suite.addTestSuite(CancerConvergenceTest.class);
    suite.addTestSuite(SomaticCallerConfigurationTest.class);
    suite.addTestSuite(CodeCrossTest.class);
    suite.addTestSuite(CombinedPriorsComplexTest.class);
    suite.addTestSuite(CombinedPriorsSnpTest.class);
    suite.addTestSuite(CombinedPriorsTest.class);
    suite.addTestSuite(HypothesesCancerTest.class);
    suite.addTestSuite(ModelCancerContaminationTest.class);
    suite.addTestSuite(ModelCancerFactoryTest.class);
    suite.addTestSuite(PosteriorContaminatedTest.class);
    suite.addTestSuite(PosteriorPureTest.class);
    suite.addTestSuite(SomaticCliTest.class);
    suite.addTestSuite(SomaticNanoTest.class);
    suite.addTestSuite(SomaticStatisticsTest.class);
    suite.addTestSuite(TopScorerTest.class);
    suite.addTestSuite(TotalScorerTest.class);

    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
