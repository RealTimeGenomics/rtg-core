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
package com.rtg.sam;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.sam");
    suite.addTestSuite(BadSuperCigarExceptionTest.class);
    suite.addTestSuite(Sam2BamTest.class);
    suite.addTestSuite(BamIndexerTest.class);
    suite.addTestSuite(BamIndexMergeTest.class);
    suite.addTestSuite(BamIndexReaderTest.class);
    suite.addTestSuite(BamReaderTest.class);
    suite.addTestSuite(BgzfInputStreamTest.class);
    suite.addTestSuite(CappedConcurrentLinkedListTest.class);
    suite.addTestSuite(CircularBufferMultifileSinglePassReaderWindowSyncTest.class);
    suite.addTestSuite(CircularBufferMultifileSinglePassReaderWindowTest.class);
    suite.addTestSuite(DedupifyingIteratorTest.class);
    suite.addTestSuite(DedupifyingRecordIteratorTest.class);
    suite.addTestSuite(DefaultSamFilterTest.class);
    suite.addTestSuite(DuplicateSamFilterTest.class);
    suite.addTestSuite(FlushLocusTest.class);
    suite.addTestSuite(PileUpTest.class);
    suite.addTestSuite(SingleMappedParamsTest.class);
    suite.addTestSuite(RandomArrayListTest.class);
    suite.addTestSuite(SamBamReaderTest.class);
    suite.addTestSuite(SamBamRecordImplTest.class);
    suite.addTestSuite(SamPicardTest.class);
    suite.addTestSuite(SamRecordExceptionTest.class);
    suite.addTestSuite(SamCommandHelperTest.class);
    suite.addTestSuite(SamClosedFileReaderTest.class);
    suite.addTestSuite(SamCompareUtilsTest.class);
    suite.addTestSuite(SamFileAndRecordTest.class);
    suite.addTestSuite(SamFilterOptionsTest.class);
    suite.addTestSuite(SamFilterParamsTest.class);
    suite.addTestSuite(SamIteratorTaskTest.class);
    suite.addTestSuite(SamRegionRestrictionTest.class);
    suite.addTestSuite(SamRestrictingIteratorTest.class);
    suite.addTestSuite(MappedParamsTest.class);
    suite.addTestSuite(MultifileIteratorTest.class);
    suite.addTestSuite(SamFilterChainTest.class);
    suite.addTestSuite(SamValidatorTest.class);
    suite.addTestSuite(SamValidatorCgHelperTest.class);
    suite.addTestSuite(SamValidatorCliTest.class);
    suite.addTestSuite(SkipInvalidRecordsIteratorTest.class);
    suite.addTestSuite(SamUtilsTest.class);
    suite.addTestSuite(SamRangeUtilsTest.class);
    suite.addTestSuite(SamMergeCliTest.class);
    suite.addTestSuite(CigarFormatterTest.class);
    suite.addTestSuite(MismatchPositionFormatterTest.class);
    suite.addTestSuite(SuperCigarTest.class);
    suite.addTestSuite(SuperCigarParserTest.class);
    suite.addTestSuite(SuperCigarValidatorTest.class);
    suite.addTestSuite(SimulatedSuperCigarUnrollerTest.class);
    suite.addTestSuite(ReadGroupUtilsTest.class);
    suite.addTestSuite(ThreadedMultifileIteratorTest.class);
    suite.addTestSuite(ThreadedMultifileIteratorWrapperTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
