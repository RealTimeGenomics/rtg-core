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
package com.rtg.util.intervals;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.util.intervals");
    suite.addTestSuite(StatusIntervalTest.class);
    suite.addTestSuite(RangeTest.class);
    suite.addTestSuite(LongRangeTest.class);
    suite.addTestSuite(SequenceIdLocusComparatorTest.class);
    suite.addTestSuite(SequenceIdLocusSimpleTest.class);
    suite.addTestSuite(SequenceNameLocusComparatorTest.class);
    suite.addTestSuite(SequenceNameLocusSimpleTest.class);
    suite.addTestSuite(RangeListTest.class);
    suite.addTestSuite(ReferenceRangesTest.class);
    suite.addTestSuite(ReferenceRegionsTest.class);
    suite.addTestSuite(RegionRestrictionTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

