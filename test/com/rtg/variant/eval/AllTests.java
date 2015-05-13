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
package com.rtg.variant.eval;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant.eval");
    suite.addTestSuite(DetectedVariantTest.class);
    suite.addTestSuite(EvalSynchronizerTest.class);
    suite.addTestSuite(HalfPathTest.class);
    suite.addTestSuite(VcfEvalTaskTest.class);
    suite.addTestSuite(VcfEvalCliTest.class);
    suite.addTestSuite(VcfEvalParamsTest.class);
    suite.addTestSuite(OrientedVariantTest.class);
    suite.addTestSuite(PathTest.class);
    suite.addTestSuite(RocContainerTest.class);
    suite.addTestSuite(RocFilterTest.class);
    suite.addTestSuite(RocLineTest.class);
    suite.addTestSuite(RocSlopeTest.class);
    suite.addTestSuite(RocScoreFieldTest.class);
    suite.addTestSuite(RocSortOrderTest.class);
    suite.addTestSuite(PhasingEvaluatorTest.class);
    suite.addTestSuite(HaplotypePlaybackTest.class);
    suite.addTestSuite(VcfRecordTabixCallableTest.class);
    suite.addTestSuite(TabixVcfRecordSetTest.class);
    suite.addTestSuite(VariantSetTypeTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
