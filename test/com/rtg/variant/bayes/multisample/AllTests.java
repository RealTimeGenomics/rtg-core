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
package com.rtg.variant.bayes.multisample;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant.bayes.multisample");

    suite.addTest(com.rtg.variant.bayes.multisample.cancer.AllTests.suite());
    suite.addTest(com.rtg.variant.bayes.multisample.family.AllTests.suite());
    suite.addTest(com.rtg.variant.bayes.multisample.forwardbackward.AllTests.suite());
    suite.addTest(com.rtg.variant.bayes.multisample.lineage.AllTests.suite());
    suite.addTest(com.rtg.variant.bayes.multisample.multithread.AllTests.suite());
    suite.addTest(com.rtg.variant.bayes.multisample.population.AllTests.suite());
    suite.addTest(com.rtg.variant.bayes.multisample.singleton.AllTests.suite());
    suite.addTest(com.rtg.variant.bayes.multisample.statistics.AllTests.suite());

    //$JUnit-BEGIN$
    suite.addTestSuite(BedComplexitiesWriterTest.class);
    suite.addTestSuite(ChunkInfoTest.class);
    suite.addTestSuite(ComplexCallerTest.class);
    suite.addTestSuite(ComplexitiesTest.class);
    suite.addTestSuite(ComplexRegionTest.class);
    suite.addTestSuite(EquivalentFilterTest.class);
    suite.addTestSuite(HypothesisScoreTest.class);
    suite.addTestSuite(HypothesisScoresTest.class);
    suite.addTestSuite(HaploidDiploidHypothesesTest.class);
    suite.addTestSuite(IndividualSampleProcessorTest.class);
    suite.addTestSuite(IndividualSampleFactoryTest.class);
    suite.addTestSuite(MultisampleTaskTest.class);
    suite.addTestSuite(MultisampleTaskSequentialTest.class);
    suite.addTestSuite(MultisampleTaskT1Test.class);
    suite.addTestSuite(MultisampleTaskT4Test.class);
    suite.addTestSuite(MultisampleUtilsTest.class);
    suite.addTestSuite(OutputUtilsTest.class);
    suite.addTestSuite(PriorContainerTest.class);
    suite.addTestSuite(RepeatAntiMeasurerTest.class);
    suite.addTestSuite(SimpleRepeatMeasurerTest.class);
    suite.addTestSuite(SingleNMerRepeatMeasurerTest.class);
    suite.addTestSuite(UtilsTest.class);
    //$JUnit-END$
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
