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
package com.rtg.segregation;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.segregation");

    suite.addTestSuite(CrossOverTest.class);
    suite.addTestSuite(FamilyGtTest.class);
    suite.addTestSuite(GTypeMendelianFactoryTest.class);
    suite.addTestSuite(GTypeTest.class);
    suite.addTestSuite(InnerAppendTest.class);
    suite.addTestSuite(LabellingTest.class);
    suite.addTestSuite(LinkedSetTest.class);
    suite.addTestSuite(MismatchingPloidyExceptionTest.class);
    suite.addTestSuite(PatternArrayTest.class);
    suite.addTestSuite(PatternDiffTest.class);
    suite.addTestSuite(PatternHolderTest.class);
    suite.addTestSuite(PatternTest.class);
    suite.addTestSuite(SearchContainerTest.class);
    suite.addTestSuite(SearchTypeTest.class);
    suite.addTestSuite(SegregationBlockTest.class);
    suite.addTestSuite(SegregationCheckerCliTest.class);
    suite.addTestSuite(SegregationVcfSearchTest.class);
    suite.addTestSuite(WhichParentTest.class);
    return suite;
  }
}
