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
package com.rtg.protein;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.protein");
    suite.addTestSuite(GotohProteinEditDistanceTest.class);
    suite.addTestSuite(MapXCliTest.class);
    suite.addTestSuite(MapXFunctionalTest.class);
    suite.addTestSuite(MapXValidatorTest.class);
    suite.addTestSuite(ProteinAlignmentResultTest.class);
    suite.addTestSuite(ProteinBitScorerTest.class);
    suite.addTestSuite(ProteinOutputProcessorTest.class);
    suite.addTestSuite(ProteinReadIndexerTest.class);
    suite.addTestSuite(ReadLengthHashingStateTest.class);
    suite.addTestSuite(ScoringHelperTest.class);
    suite.addTestSuite(SequenceLengthBucketsTest.class);
    suite.addTestSuite(SharedProteinResourcesTest.class);
    suite.addTestSuite(SharedStatusCollectorTest.class);
    suite.addTestSuite(TopEqualProteinOutputProcessorTest.class);
    suite.addTestSuite(TopEqualProteinImplementationTest.class);
    suite.addTestSuite(TopNProteinImplementationTest.class);
    suite.addTestSuite(TopNProteinOutputProcessorTest.class);
    suite.addTestSuite(MapxRenameTest.class);
    suite.addTestSuite(MapXStatisticsTest.class);
    return suite;
  }
}
