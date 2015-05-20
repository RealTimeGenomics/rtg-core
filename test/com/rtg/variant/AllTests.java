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
package com.rtg.variant;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant");
    suite.addTestSuite(BaseQualityPhredScalerTest.class);
    suite.addTestSuite(BaseQualityReadPositionPhredScalerTest.class);
    suite.addTestSuite(CalibratedPerSequenceThresholdTest.class);
    suite.addTestSuite(CalibratedMachineErrorChooserTest.class);
    suite.addTestSuite(CalibratedMachineErrorParamsTest.class);
    suite.addTestSuite(DefaultMachineErrorChooserTest.class);
    suite.addTestSuite(GenomeConnectivityTest.class);
    suite.addTestSuite(GenomePriorParamsTest.class);
    suite.addTestSuite(MachineErrorParamsTest.class);
    suite.addTestSuite(NoQualityExceptionTest.class);
    suite.addTestSuite(PosteriorUtilsTest.class);
    suite.addTestSuite(ReadPositionPhredScalerTest.class);
    suite.addTestSuite(ReadGroupMachineErrorChooserTest.class);
    suite.addTestSuite(RecalibratingSamRecordPopulatorTest.class);
    suite.addTestSuite(SamRecordPopulatorTest.class);
    suite.addTestSuite(SamToMatchCigarTest.class);
    suite.addTestSuite(StaticThresholdTest.class);
    suite.addTestSuite(ThreadingEnvironmentTest.class);
    suite.addTestSuite(VariantTest.class);
    suite.addTestSuite(VariantAlignmentRecordPopulatorTest.class);
    suite.addTestSuite(VariantAlignmentRecordTest.class);
    suite.addTestSuite(VariantOutputLevelTest.class);
    suite.addTestSuite(VariantParamsTest.class);
    suite.addTestSuite(VariantStatisticsTest.class);
    suite.addTestSuite(VariantNanoTest.class);
    suite.addTestSuite(VariantLocusTest.class);
    suite.addTestSuite(VariantSampleTest.class);
    suite.addTestSuite(VariantSampleInfoTest.class);
    suite.addTestSuite(VcfStatsCliTest.class);

    suite.addTest(com.rtg.variant.bayes.AllTests.suite());
    suite.addTest(com.rtg.variant.cnv.AllTests.suite());
    suite.addTest(com.rtg.variant.coverage.AllTests.suite());
    suite.addTest(com.rtg.variant.dna.AllTests.suite());
    suite.addTest(com.rtg.variant.match.AllTests.suite());
    suite.addTest(com.rtg.variant.realign.AllTests.suite());
    suite.addTest(com.rtg.variant.eval.AllTests.suite());
    suite.addTest(com.rtg.variant.sv.AllTests.suite());
    suite.addTest(com.rtg.variant.util.AllTests.suite());
    suite.addTest(com.rtg.variant.format.AllTests.suite());
    suite.addTest(com.rtg.variant.avr.AllTests.suite());

    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
