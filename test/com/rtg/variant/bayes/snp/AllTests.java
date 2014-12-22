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
package com.rtg.variant.bayes.snp;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant.bayes.snp");

    suite.addTestSuite(CigarParserModelTest.class);
    suite.addTestSuite(DescriptionCommonTest.class);
    suite.addTestSuite(DescriptionNoneTest.class);
    suite.addTestSuite(DescriptionSnpTest.class);
    suite.addTestSuite(EvidenceIndelFactoryTest.class);
    suite.addTestSuite(EvidenceIndelTest.class);
    suite.addTestSuite(EvidenceQTest.class);
    suite.addTestSuite(EvidenceQFactoryTest.class);
    suite.addTestSuite(HypothesesNoneTest.class);
    suite.addTestSuite(HypothesesSnpTest.class);
    suite.addTestSuite(HypothesesPriorTest.class);
    suite.addTestSuite(ReferenceBasedBufferTest.class);
    suite.addTestSuite(SwitchingReferenceBasedBufferTest.class);
    suite.addTestSuite(IndelDetectorFactoryTest.class);
    suite.addTestSuite(IndelDetectorTest.class);
    suite.addTestSuite(EvidenceMatcherTest.class);
    suite.addTestSuite(IndelMatcherTest.class);
    suite.addTestSuite(ModelNoneFactoryTest.class);
    suite.addTestSuite(ModelSnpFactoryTest.class);
    suite.addTestSuite(StatisticsSnpTest.class);

    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
