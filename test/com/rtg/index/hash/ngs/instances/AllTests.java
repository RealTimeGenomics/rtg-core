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
package com.rtg.index.hash.ngs.instances;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory.
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.index.hash.ngs.instances");

    suite.addTestSuite(DummyCGMaskTest.class);
    suite.addTestSuite(CGMaska0Test.class);
    suite.addTestSuite(CGMaska1b1Test.class);
    suite.addTestSuite(CGMaska1b1altTest.class);
    suite.addTestSuite(CGMaska15b1Test.class);
    suite.addTestSuite(CGMaska15b1altTest.class);

    suite.addTestSuite(MaskL51w12s3e3Test.class);
    suite.addTestSuite(MaskL51w15s3e2Test.class);
    suite.addTestSuite(MaskL51w17s2e2Test.class);
    suite.addTestSuite(MaskL51w18s2e1Test.class);
    suite.addTestSuite(MaskL36w18s3e1Test.class);
    suite.addTestSuite(MaskL35w15s2e1Test.class);
    suite.addTestSuite(MaskL30w12s3e1Test.class);

    suite.addTestSuite(SplitL36w12s2e2Test.class);
    suite.addTestSuite(SplitL36w18s1e1Test.class);
    suite.addTestSuite(SplitL36w18s2e1Test.class);
    suite.addTestSuite(SplitL4w2s1e1bTest.class);
    suite.addTestSuite(SplitL4w2s1e1Test.class);
    suite.addTestSuite(SplitL4w4s0e0Test.class);
    suite.addTestSuite(SplitL36w12s3e3Test.class);
    suite.addTestSuite(SplitL36w12s3e2Test.class);

    suite.addTestSuite(SubstituteTest.class);
    suite.addTestSuite(SubstituteIndelTest.class);

    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

