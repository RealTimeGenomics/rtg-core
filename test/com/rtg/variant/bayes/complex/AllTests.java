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
package com.rtg.variant.bayes.complex;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant.bayes.snp");
    suite.addTestSuite(ComplexTemplateTest.class);
    suite.addTestSuite(DescriptionComplexTest.class);
    suite.addTestSuite(EvidenceComplexTest.class);
    suite.addTestSuite(HypothesesComplexTest.class);
    suite.addTestSuite(IndelUtilsTest.class);
    suite.addTestSuite(InsertionPriorTest.class);
    suite.addTestSuite(IonTorrentCallFilterTest.class);
    suite.addTestSuite(KappaImplementationTest.class);
    suite.addTestSuite(KappaMemoTest.class);
    suite.addTestSuite(LineageDenovoCheckerTest.class);
    suite.addTestSuite(MatchMultiSetTest.class);
    suite.addTestSuite(MaxShiftCigarParserTest.class);
    suite.addTestSuite(MendelianDenovoCheckerTest.class);
    suite.addTestSuite(ScoreInterfaceMemoTest.class);
    suite.addTestSuite(ScoreInterfaceMemoIonTTest.class);
    suite.addTestSuite(SingleCountsTest.class);
    suite.addTestSuite(StatisticsComplexTest.class);
    suite.addTestSuite(TauImplementationTest.class);
    suite.addTestSuite(TrimmingTest.class);
    return suite;
  }
}
