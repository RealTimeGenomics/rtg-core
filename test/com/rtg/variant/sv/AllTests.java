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
package com.rtg.variant.sv;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant.sv");
    suite.addTestSuite(AllCountsTest.class);
    suite.addTestSuite(BayesianSignalTest.class);
    suite.addTestSuite(BreakpointBayesianSignalTest.class);
    suite.addTestSuite(CorrectionsTest.class);
    suite.addTestSuite(CumulativeSamCountsTest.class);
    suite.addTestSuite(DeleteBoundaryBayesianSignalTest.class);
    suite.addTestSuite(DistributionArrayTest.class);
    suite.addTestSuite(DistributionConstantTest.class);
    suite.addTestSuite(DistributionStepTest.class);
    suite.addTestSuite(DistributionUtilsTest.class);
    suite.addTestSuite(DuplicateDonorBayesianSignalTest.class);
    suite.addTestSuite(HeterozygousBayesianSignalTest.class);
    suite.addTestSuite(NormalBayesianSignalTest.class);
    suite.addTestSuite(NovelInsertionBayesianSignalTest.class);
    suite.addTestSuite(ReadGroupStateTest.class);
    suite.addTestSuite(ReadGroupStateTest.class);
    suite.addTestSuite(ReadGroupStatsTest.class);
    suite.addTestSuite(ReadGroupStatsCalculatorTest.class);
    suite.addTestSuite(SamArrayTest.class);
    suite.addTestSuite(SignalCountTest.class);
    suite.addTestSuite(SignalConstantLnTest.class);
    suite.addTestSuite(SignalDistributionLnTest.class);
    suite.addTestSuite(SignalSumTest.class);
    suite.addTestSuite(SvCliUtilsTest.class);
    suite.addTestSuite(SvInterestingRegionExtractorTest.class);
    suite.addTestSuite(SvParamsTest.class);
    suite.addTestSuite(SvToolCliTest.class);
    suite.addTestSuite(SvToolTaskTest.class);
    suite.addTestSuite(SvToolParamsTest.class);
    suite.addTestSuite(UnmatedAugmenterTest.class);
    suite.addTestSuite(UnmatedAugmenterCliTest.class);

    suite.addTest(com.rtg.variant.sv.discord.AllTests.suite());
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
