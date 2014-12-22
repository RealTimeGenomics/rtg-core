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
package com.rtg.vcf;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.vcf");

    //$JUnit-BEGIN$
    suite.addTest(com.rtg.vcf.annotation.AllTests.suite());
    suite.addTest(com.rtg.vcf.header.AllTests.suite());
    suite.addTest(com.rtg.vcf.validator.AllTests.suite());
    suite.addTestSuite(AdjacencyTest.class);
    suite.addTestSuite(ChildPhasingVcfAnnotatorTest.class);
    suite.addTestSuite(SegregationVcfAnnotatorTest.class);
    suite.addTestSuite(VcfFilterStatisticsTest.class);
    suite.addTestSuite(SnpIntersectionTest.class);
    suite.addTestSuite(VcfAnnotatorCliTest.class);
    suite.addTestSuite(VcfFilterCliTest.class);
    suite.addTestSuite(VcfFilterStripperTest.class);
    suite.addTestSuite(VcfFormatDoubleAnnotatorTest.class);
    suite.addTestSuite(VcfFormatStringAnnotatorTest.class);
    suite.addTestSuite(VcfFormatStripperTest.class);
    suite.addTestSuite(VcfInfoDoubleAnnotatorTest.class);
    suite.addTestSuite(VcfInfoFilterTest.class);
    suite.addTestSuite(VcfInfoIntegerAnnotatorTest.class);
    suite.addTestSuite(VcfInfoStripperTest.class);
    suite.addTestSuite(VcfInfoPerAltIntegerAnnotatorTest.class);
    suite.addTestSuite(VcfMergeTest.class);
    suite.addTestSuite(VcfReaderTest.class);
    suite.addTestSuite(VcfRecordTest.class);
    suite.addTestSuite(VcfSampleStripperTest.class);
    suite.addTestSuite(VcfSubsetTest.class);
    suite.addTestSuite(VcfUtilsTest.class);
    suite.addTestSuite(VcfWriterTest.class);
    //$JUnit-END$

    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
