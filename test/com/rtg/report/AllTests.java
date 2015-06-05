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

package com.rtg.report;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite(AllTests.class.getName());
    suite.addTestSuite(CombinedReportTest.class);
    suite.addTestSuite(MapFSummaryReportTest.class);
    suite.addTestSuite(MapReportTest.class);
    suite.addTestSuite(MapReportTemplateTest.class);
    suite.addTestSuite(MapXSummaryReportTest.class);
    suite.addTestSuite(MetagenomicsReportTemplateTest.class);
    suite.addTestSuite(ReportableModuleTest.class);
    suite.addTestSuite(ReportCliTest.class);
    suite.addTestSuite(ReportTemplateTest.class);
    suite.addTestSuite(ReportTypeTest.class);
    suite.addTestSuite(RtgVelocityLogChuteTest.class);
    suite.addTestSuite(SimilarityReportTest.class);
    suite.addTestSuite(SpeciesReportTest.class);
    suite.addTestSuite(VelocityReportUtilsTest.class);
    return suite;
  }

}
