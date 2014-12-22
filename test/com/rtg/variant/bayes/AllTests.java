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
package com.rtg.variant.bayes;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant.bayes");

    suite.addTest(com.rtg.variant.bayes.complex.AllTests.suite());
    suite.addTest(com.rtg.variant.bayes.multisample.AllTests.suite());
    suite.addTest(com.rtg.variant.bayes.snp.AllTests.suite());

    suite.addTestSuite(ArrayGenotypeMeasureTest.class);
    suite.addTestSuite(BinaryFactorTest.class);
    suite.addTestSuite(CodeDiploidTest.class);
    suite.addTestSuite(CodeHaploidTest.class);
    suite.addTestSuite(AlleleStatisticsDoubleTest.class);
    suite.addTestSuite(AlleleStatisticsIntTest.class);
    suite.addTestSuite(DefaultSortedFactorTest.class);
    suite.addTestSuite(EntropyBinomialTest.class);
    suite.addTestSuite(EntropyRandomTest.class);
    suite.addTestSuite(EvidenceTest.class);
    suite.addTestSuite(NormalizedFactorTest.class);
    suite.addTestSuite(HypothesesTest.class);
    suite.addTestSuite(ModelNoneTest.class);
    suite.addTestSuite(ModelTest.class);
    suite.addTestSuite(NoNonIdentityMeasureTest.class);
    suite.addTestSuite(StatisticsIntTest.class);
    suite.addTestSuite(StatisticsDoubleTest.class);
    suite.addTestSuite(UnitFactorTest.class);

    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
