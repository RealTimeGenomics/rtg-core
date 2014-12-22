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
package com.rtg.index;


import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory. Run from the command
 * line with:<p>
 *
 * java -ea com.reeltwo.cartesian.AllTests
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.index");

    suite.addTest(com.rtg.index.hash.AllTests.suite());
    suite.addTest(com.rtg.index.params.AllTests.suite());
    suite.addTest(com.rtg.index.queue.AllTests.suite());
    suite.addTest(com.rtg.index.similarity.AllTests.suite());

    suite.addTestSuite(SearchUtilsTest.class);
    suite.addTestSuite(BigBitVectorTest.class);
    suite.addTestSuite(BitVectorTest.class);
    suite.addTestSuite(FinderTest.class);
    suite.addTestSuite(HashBitHandleTest.class);
    suite.addTestSuite(HashBitVectorTest.class);
    suite.addTestSuite(IndexBaseTest.class);
    suite.addTestSuite(IndexCompressedExtendedTest.class);
    suite.addTestSuite(IndexCompressedTest.class);
    suite.addTestSuite(IndexSimpleTest.class);
    suite.addTestSuite(IntSetTest.class);
    suite.addTestSuite(IntSetSingleTest.class);
    suite.addTestSuite(IntSetWindowTest.class);
    suite.addTestSuite(IndexSetTest.class);
    suite.addTestSuite(IndexSorterTest.class);
    suite.addTestSuite(IndexUtilsTest.class);
    suite.addTestSuite(ReadIndexCoverageUtilsTest.class);
    suite.addTestSuite(SparseFrequencyHistogramTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

