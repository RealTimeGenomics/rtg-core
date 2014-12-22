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
package com.rtg.visualization;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.visualization");
    suite.addTestSuite(AviewTest.class);
    suite.addTestSuite(AviewModelTest.class);
    suite.addTestSuite(AviewParamsTest.class);
    suite.addTestSuite(AviewParamsBuilderTest.class);
    suite.addTestSuite(AviewVariantTest.class);
    suite.addTestSuite(ReferenceHelperTest.class);
    suite.addTestSuite(CigarHelperTest.class);
    suite.addTestSuite(DisplayHelperTest.class);
    suite.addTestSuite(AnsiDisplayHelperTest.class);
    suite.addTestSuite(HtmlDisplayHelperTest.class);
    suite.addTestSuite(SamAssistanceCgTest.class);
    suite.addTestSuite(SamAssistanceCgLegacyTest.class);
    suite.addTestSuite(SamAssistanceSimpleTest.class);
    suite.addTestSuite(SamHelperTest.class);
    suite.addTestSuite(SuperCigarUnrollerTest.class);
    suite.addTestSuite(UnrolledReadTest.class);
    suite.addTestSuite(VariantHelperTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
