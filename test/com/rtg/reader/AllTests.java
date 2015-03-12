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
package com.rtg.reader;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * AllTests for reader package
 */
public class AllTests extends TestCase {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.reader");
    suite.addTest(CompressedMemorySequencesReaderTest.suite());
    suite.addTest(DefaultSequencesReaderTest.suite());
    suite.addTest(IndexFileTest.suite());
    suite.addTest(FormatCliTest.suite());
    suite.addTest(ReverseComplementingReaderTest.suite());
    suite.addTest(SequenceStreamManagerTest.suite());
    suite.addTest(SequencesWriterTest.suite());

    suite.addTestSuite(AlternatingSequencesReaderTest.class);
    suite.addTestSuite(AlternatingSequencesWriterTest.class);
    suite.addTestSuite(ArrayPrereadNamesTest.class);
    suite.addTestSuite(BestSumReadTrimmerTest.class);
    suite.addTestSuite(SdfFileUtilsTest.class);
    suite.addTestSuite(Cg2SdfTest.class);
    suite.addTestSuite(CompressedMemorySequencesReader2Test.class);
    suite.addTestSuite(CompressedMemorySequencesWriterTest.class);
    suite.addTestSuite(ConcatSequenceDataSourceTest.class);
    suite.addTestSuite(CorruptSdfExceptionTest.class);
    suite.addTestSuite(DataFileOpenerFactoryTest.class);
    suite.addTestSuite(DataInMemoryTest.class);
    suite.addTestSuite(EmptyStringPrereadNamesTest.class);
    suite.addTestSuite(FastaUtilsTest.class);
    suite.addTestSuite(FastaSequenceDataSourceTest.class);
    suite.addTestSuite(FastqSequenceDataSourceTest.class);
    suite.addTestSuite(FileBitwiseInputStreamTest.class);
    suite.addTestSuite(FileBitwiseOutputStreamTest.class);
    suite.addTestSuite(FileCompressedInputStreamTest.class);
    suite.addTestSuite(FileCompressedOutputStreamTest.class);
    suite.addTestSuite(FileStreamIteratorTest.class);
    suite.addTestSuite(GilletteReadTrimmerTest.class);
    suite.addTestSuite(InputFormatTest.class);
    suite.addTestSuite(LabelTest.class);
    suite.addTestSuite(MappedSamBamSequenceDataSourceTest.class);
    suite.addTestSuite(NameDuplicateDetectorTest.class);
    suite.addTestSuite(NameFilePairTest.class);
    suite.addTestSuite(PointerFileHandlerTest.class);
    suite.addTestSuite(PointerFileLookupTest.class);
    suite.addTestSuite(PrereadArmTest.class);
    suite.addTestSuite(PrereadHashFunctionTest.class);
    suite.addTestSuite(PrereadNamesTest.class);
    suite.addTestSuite(SdfSplitterTest.class);
    suite.addTestSuite(SdfStatisticsTest.class);
    suite.addTestSuite(Sdf2FastaTest.class);
    suite.addTestSuite(Sdf2FastqTest.class);
    suite.addTestSuite(Sdf2QualaTest.class);
    suite.addTestSuite(PrereadTypeTest.class);
    suite.addTestSuite(PrereadVerifierTest.class);
    suite.addTestSuite(ReadHelperTest.class);
    suite.addTestSuite(ReaderUtilsTest.class);
    suite.addTestSuite(RightSimplePrereadNamesTest.class);
    suite.addTestSuite(SamBamSequenceDataSourceTest.class);
    suite.addTestSuite(SdfFilterTest.class);
    suite.addTestSuite(SdfIdTest.class);
    suite.addTestSuite(SdfSubseqTest.class);
    suite.addTestSuite(SdfSubsetTest.class);
    suite.addTestSuite(SdfUtilsTest.class);
    suite.addTestSuite(SequenceDataLoaderTest.class);
    suite.addTestSuite(SequencesReaderFactoryTest.class);
    suite.addTestSuite(SimplePrereadNamesTest.class);
    suite.addTestSuite(SourceTemplateReadWriterTest.class);
    suite.addTestSuite(TsvSequenceDataSourceTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

