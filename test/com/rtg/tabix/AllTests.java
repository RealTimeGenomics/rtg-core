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
package com.rtg.tabix;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.sam");
    suite.addTestSuite(AlleleCountsPositionReaderTest.class);
    suite.addTestSuite(BgZipTest.class);
    suite.addTestSuite(BlockCompressedLineReaderTest.class);
    suite.addTestSuite(BrLineReaderTest.class);
    suite.addTestSuite(ExtractCliTest.class);
    suite.addTestSuite(GenericPositionReaderTest.class);
    suite.addTestSuite(IndexerCliTest.class);
    suite.addTestSuite(IndexingStreamCreatorTest.class);
    suite.addTestSuite(IndexUtilsTest.class);
    suite.addTestSuite(TabixIndexMergeTest.class);
    suite.addTestSuite(SamPositionReaderTest.class);
    suite.addTestSuite(SequenceIndexContainerTest.class);
    suite.addTestSuite(TabixLineReaderTest.class);
    suite.addTestSuite(TabixHeaderTest.class);
    suite.addTestSuite(TabixIndexerTest.class);
    suite.addTestSuite(TabixIndexReaderTest.class);
    suite.addTestSuite(VcfPositionReaderTest.class);

    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
