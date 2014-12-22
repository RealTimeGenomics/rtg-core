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
package com.rtg.position.output;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory. Run from the command
 * line with:<p>
 *
 * java com.reeltwo.cartesian.AllTests
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.position.output");

    suite.addTestSuite(Blosum62ScoresTest.class);
    suite.addTestSuite(Blosum62Test.class);
    suite.addTestSuite(DummyPositionOutputTest.class);
    suite.addTestSuite(GapBucketsTest.class);
    suite.addTestSuite(GapBucketsInfoTest.class);
    suite.addTestSuite(GappedDistributionTest.class);
    suite.addTestSuite(GappedOutputTest.class);
    suite.addTestSuite(GapProbabilitiesScorerTest.class);
    suite.addTestSuite(GappedRegionTest.class);
    suite.addTestSuite(MismatchScoresTest.class);
    suite.addTestSuite(NgsWordOutputTest.class);
    suite.addTestSuite(OffsetTest.class);
    suite.addTestSuite(OutputFormatTypeTest.class);
    suite.addTestSuite(OutputUtilsTest.class);
    suite.addTestSuite(PositionDistributionParamsTest.class);
    suite.addTestSuite(PositionOutputTest.class);
    suite.addTestSuite(PositionOutputParamsTest.class);
    suite.addTestSuite(PositionParamsTest.class);
    suite.addTestSuite(RawPositionWriterTest.class);
    suite.addTestSuite(SegmentTest.class);
    suite.addTestSuite(SegmentCollectionTest.class);
    suite.addTestSuite(SegmentCollectorTest.class);
    suite.addTestSuite(SegmentOutputTest.class);
    suite.addTestSuite(ScannerTest.class);
    suite.addTestSuite(SurrogateGappedRegionTest.class);
    suite.addTestSuite(UtilityTest.class);
    suite.addTestSuite(WordOutputTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
