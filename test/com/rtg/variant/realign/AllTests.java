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
package com.rtg.variant.realign;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant.realign");
    suite.addTestSuite(ByteArrayAdaptorTest.class);
    suite.addTestSuite(CgRealignTest.class);
    suite.addTestSuite(EnvironmentCombinedTest.class);
    suite.addTestSuite(EnvironmentDecoratorTest.class);
    suite.addTestSuite(EnvironmentHomopolymerTest.class);
    suite.addTestSuite(EnvironmentImplementationTest.class);
    suite.addTestSuite(EnvironmentSNPTest.class);
    suite.addTestSuite(EnvironmentSubreadTest.class);
    suite.addTestSuite(HomopolymerMatrixTest.class);
    suite.addTestSuite(HomoPolymerParamsTest.class);
    suite.addTestSuite(HomopolymerRepeatsTest.class);
    suite.addTestSuite(InvertedEnvironmentTest.class);
    suite.addTestSuite(InvertCgTemplateEnvironmentTest.class);
    suite.addTestSuite(DeltaImplementationTest.class);
    suite.addTestSuite(RealignParamsImplementationTest.class);
    suite.addTestSuite(DeltaSlowlyTest.class);
    suite.addTestSuite(RealignParamsGenomeTest.class);
    suite.addTestSuite(AlignmentEnvironmentReadTest.class);
    suite.addTestSuite(AlignmentEnvironmentGenomeTest.class);
    suite.addTestSuite(AlignmentEnvironmentGenomeSubstitutionTest.class);
    suite.addTestSuite(AlignmentEnvironmentCGTest.class);
    suite.addTest(ScoreFastUnderflowCGTest.suite());
    suite.addTestSuite(ScoreFastUnderflowHomopolymerTest.class);
    suite.addTest(ScoreFastUnderflowTest.suite());
    suite.addTest(ScoreMatrixCGReverseTest.suite());
    suite.addTest(ScoreMatrixCGTest.suite());
    suite.addTest(ScoreMatrixReverseTest.suite());
    suite.addTest(ScoreMatrixTest.suite());

    suite.addTest(com.rtg.variant.realign.treealign.AllTests.suite());

    return suite;
  }
}
