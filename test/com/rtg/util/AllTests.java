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
package com.rtg.util;


import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory.
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.util");
    suite.addTest(com.rtg.util.arithcode.AllTests.suite());
    suite.addTest(com.rtg.util.array.AllTests.suite());
    suite.addTest(com.rtg.util.cli.AllTests.suite());
    suite.addTest(com.rtg.util.diagnostic.AllTests.suite());
    suite.addTest(com.rtg.util.gzip.AllTests.suite());
    suite.addTest(com.rtg.util.io.AllTests.suite());
    suite.addTest(com.rtg.util.format.AllTests.suite());
    suite.addTest(com.rtg.util.integrity.AllTests.suite());
    suite.addTest(com.rtg.util.intervals.AllTests.suite());
    suite.addTest(com.rtg.util.iterators.AllTests.suite());
    suite.addTest(com.rtg.util.machine.AllTests.suite());
    suite.addTest(com.rtg.util.memo.AllTests.suite());
    suite.addTest(com.rtg.util.memory.AllTests.suite());
    suite.addTest(com.rtg.util.store.AllTests.suite());
    suite.addTest(com.rtg.util.test.AllTests.suite());
    suite.addTest(com.rtg.util.bytecompression.AllTests.suite());
    suite.addTest(EnvironmentTest.suite());
    suite.addTest(PairTest.suite());
    suite.addTest(SpawnJvmTest.suite());

    suite.addTestSuite(AutoAddMapTest.class);
    suite.addTestSuite(BitPack2IntoLongTest.class);
    suite.addTestSuite(BitPack3IntoIntTest.class);
    suite.addTestSuite(BitPack3IntoLongTest.class);
    suite.addTestSuite(BitPackHelperIntTest.class);
    suite.addTestSuite(BitPackHelperLongTest.class);
    suite.addTestSuite(BoundedDoubleTest.class);
    suite.addTestSuite(ByteUtilsTest.class);
    suite.addTestSuite(ChiSquaredTest.class);
    suite.addTestSuite(ChooseMemoryTest.class);
    suite.addTestSuite(ComparableComparatorTest.class);
    suite.addTestSuite(ComparablePairTest.class);
    suite.addTestSuite(CompareHelperTest.class);
    suite.addTestSuite(ConstantsTest.class);
    suite.addTestSuite(DummyEncodingTest.class);
    suite.addTestSuite(DummyUnauthorizedAccessExceptionTest.class);
    suite.addTestSuite(EnumHelperTest.class);
    suite.addTestSuite(EthernetAddressTest.class);
    suite.addTestSuite(FlexArrayTest.class);
    suite.addTestSuite(HistogramTest.class);
    suite.addTestSuite(HtmlReportHelperTest.class);
    suite.addTestSuite(IntegerOrPercentageTest.class);
    suite.addTestSuite(InvalidParamsExceptionTest.class);
    suite.addTestSuite(IORunnableProxyTest.class);
    suite.addTestSuite(LongUtilsTest.class);
    suite.addTestSuite(BasicLinkedListNodeTest.class);
    suite.addTestSuite(MapSetTest.class);
    suite.addTestSuite(MathUtilsTest.class);
    suite.addTestSuite(MaxShiftFactorTest.class);
    suite.addTestSuite(MaxShiftUtilsTest.class);
    suite.addTestSuite(MD5UtilsTest.class);
    suite.addTestSuite(MultiMapTest.class);
    suite.addTestSuite(MultiSetTest.class);
    suite.addTestSuite(DoubleMultiSetTest.class);
    suite.addTestSuite(SortedMultiSetTest.class);
    suite.addTestSuite(NormalDistributionTest.class);
    suite.addTestSuite(NullStreamUtilsTest.class);
    suite.addTestSuite(ObjectParamsTest.class);
    suite.addTestSuite(ParamsBuilderTest.class);
    suite.addTestSuite(PortableRandomTest.class);
    suite.addTestSuite(PosteriorUtilsTest.class);
    suite.addTestSuite(ProgramStateTest.class);
    suite.addTestSuite(PropertiesUtilsTest.class);
    suite.addTestSuite(QuickSortTest.class);
    suite.addTestSuite(QuickSortDoubleIntProxyTest.class);
    suite.addTestSuite(QuickSortIntIntProxyTest.class);
    suite.addTestSuite(QuickSortLongLongProxyTest.class);
    suite.addTestSuite(ReorderingQueueTest.class);
    suite.addTestSuite(ResourcesTest.class);
    suite.addTestSuite(ResultStreamHandlerTest.class);
    suite.addTestSuite(ReverseDoubleTest.class);
    suite.addTestSuite(SeparateClassLoaderTest.class);
    suite.addTestSuite(SimpleThreadPoolTest.class);
    suite.addTestSuite(SingletonPopulatorFactoryTest.class);
    suite.addTestSuite(SizeSplitTest.class);
    suite.addTestSuite(StandardDeviationTest.class);
    suite.addTestSuite(StatisticsTest.class);
    suite.addTestSuite(StringUtilsTest.class);
    suite.addTestSuite(SynchronizedLinkedListTest.class);
    suite.addTestSuite(TextTableTest.class);
    suite.addTestSuite(TsvUtilsTest.class);
    suite.addTestSuite(UtilsTest.class);
    suite.addTestSuite(WorkerThreadTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

