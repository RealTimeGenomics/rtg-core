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
package com.rtg.ngs;

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
    final TestSuite suite = new TestSuite("com.rtg.ngs");
    suite.addTestSuite(DeduplicatingNStoreTest.class);
    suite.addTestSuite(DummyMapOutputProcessorTest.class);
    suite.addTestSuite(FilterConcatIntermediateFilesTest.class);
    suite.addTestSuite(LongReadMaskParamsTest.class);
    suite.addTestSuite(MatchResultTest.class);
    suite.addTestSuite(NgsTaskTest.class);
    suite.addTestSuite(NgsTaskFunctionalTest.class);
    suite.addTestSuite(NgsThreadTest.class);
    suite.addTestSuite(NgsLongTest.class);
    suite.addTestSuite(NgsLongThreadsTest.class);
    suite.addTestSuite(OutputFilterTest.class);
    suite.addTestSuite(NgsCgTest.class);
    suite.addTestSuite(NgsParamsTest.class);
    suite.addTestSuite(NgsOutputParamsTest.class);
    suite.addTestSuite(NgsFilterParamsTest.class);
    suite.addTestSuite(NgsMaskParamsGeneralTest.class);
    suite.addTestSuite(NgsMaskParamsExplicitTest.class);
    suite.addTestSuite(NgsMaskParamsProteinTest.class);
    suite.addTestSuite(NgsPairedEndTest.class);
    suite.addTestSuite(DefaultOutputProcessorTest.class);
    suite.addTestSuite(ClippedOutputProcessorTest.class);
    suite.addTestSuite(DefaultOutputProcessorSynchTest.class);
    suite.addTestSuite(MapReportDataTest.class);
    suite.addTestSuite(PairedEndOutputProcessorTest.class);
    suite.addTestSuite(TopNPairedEndOutputProcessorSyncTest.class);

    suite.addTestSuite(TopNImplementationTest.class);
    suite.addTestSuite(UptoNStoreSyncTest.class);
    suite.addTestSuite(CgMapCliTest.class);
    suite.addTestSuite(CgMapFlagsValidatorTest.class);
    suite.addTestSuite(MapCliTest.class);
    suite.addTestSuite(RamMapFunctionalTest.class);
    suite.addTestSuite(MatedSamResultsFilterTest.class);
    suite.addTestSuite(UnmatedSamResultsFilterTest.class);
    suite.addTestSuite(SamResultsSansFilterTest.class);
    suite.addTestSuite(SamSequenceReaderParamsTest.class);
    suite.addTestSuite(SingleEndSamResultsFilterTest.class);
    suite.addTestSuite(DummySamResultsFilterTest.class);
    suite.addTestSuite(NgsResultStreamHandlerTest.class);
    suite.addTestSuite(ReadStatusTrackerTest.class);
    suite.addTestSuite(ReadStatusTrackerSyncTest.class);
    suite.addTestSuite(DummyAlignmentWriterThreadTest.class);
    suite.addTestSuite(SamRenameTest.class);
    suite.addTestSuite(SamSingleEndOutputProcessorTest.class);
    suite.addTestSuite(UnfilteredSingleEndOutputProcessorTest.class);
    suite.addTestSuite(UnfilteredPairedEndOutputProcessorTest.class);
    suite.addTestSuite(UnmappedSamAlignmentWriterTest.class);

    suite.addTestSuite(SharedResourcesTest.class);
    suite.addTestSuite(SingleEndMapStatisticsTest.class);
    suite.addTestSuite(PairedEndMapStatisticsTest.class);
    suite.addTestSuite(MapStatisticsFieldTest.class);
    suite.addTestSuite(MapStatisticsArmTest.class);
    suite.addTestSuite(NullOutputProcessorTest.class);

    suite.addTestSuite(MapFCliTest.class);
    suite.addTestSuite(MapFilterPairedMapStatisticsTest.class);
    suite.addTestSuite(MapFilterSingleEndMapStatisticsTest.class);
    suite.addTestSuite(DummySdfOutputProcessorTest.class);

    suite.addTestSuite(MapFlagsTest.class);
    suite.addTestSuite(MapParamsHelperTest.class);
    suite.addTestSuite(PairedTopRandomImplementationTest.class);
    suite.addTestSuite(PairedTopRandomImplementationSyncTest.class);
    suite.addTestSuite(SingleEndTopRandomImplementationTest.class);
    suite.addTestSuite(SingleEndTopRandomImplementationSyncTest.class);

    suite.addTest(com.rtg.ngs.blocking.AllTests.suite());
    suite.addTest(com.rtg.ngs.longread.AllTests.suite());
    suite.addTest(com.rtg.ngs.tempstage.AllTests.suite());
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

