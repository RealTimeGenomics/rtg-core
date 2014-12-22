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
package com.rtg.assembler;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Tests.
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.assembler");
    suite.addTestSuite(AlignmentIteratorTest.class);
    suite.addTestSuite(AlignmentSectionTest.class);
    suite.addTestSuite(AlignmentChainTest.class);
    suite.addTestSuite(ApplyReferenceTest.class);
    suite.addTestSuite(AsyncReadSourceTest.class);
    suite.addTestSuite(AsyncReadPoolTest.class);
    suite.addTestSuite(AssembleCliTest.class);
    suite.addTestSuite(AssembleParamsTest.class);
    suite.addTestSuite(AssembleTaskTest.class);
    suite.addTestSuite(BubbleExplorerTest.class);
    suite.addTestSuite(ByteKmerTest.class);
    suite.addTestSuite(CountAssemblyEndsTest.class);
    suite.addTestSuite(ConsensusTest.class);
    suite.addTestSuite(ConstraintCacheTest.class);
    suite.addTestSuite(ConstraintCollectorTest.class);
    suite.addTestSuite(ContigConcatenatedTest.class);
    suite.addTestSuite(ContigCollectorTest.class);
    suite.addTestSuite(ContigPositionTest.class);
    suite.addTestSuite(CorrectingMutatorTest.class);
    suite.addTestSuite(CorrectReadsCliTest.class);
    suite.addTestSuite(CorrectReadsTest.class);
    suite.addTestSuite(DeBruijnAssemblerTaskTest.class);
    suite.addTestSuite(DeBruijnAssemblerCliTest.class);
    suite.addTestSuite(DeBruijnGraphBuilderTest.class);
    suite.addTestSuite(DeBruijnParamsTest.class);
    suite.addTestSuite(ExtractPathTest.class);
    suite.addTestSuite(FilterPathsTest.class);
    suite.addTestSuite(GraphIndexTest.class);
    suite.addTestSuite(GraphAlignerTest.class);
    suite.addTestSuite(GraphCleanupTest.class);
    suite.addTestSuite(GraphFlowAlignerTest.class);
    suite.addTestSuite(GraphAlignmentTest.class);
    suite.addTestSuite(GraphMapCliTest.class);
    suite.addTestSuite(GraphMapStatisticsTest.class);
    suite.addTestSuite(GraphMapTest.class);
    suite.addTestSuite(GraphMapTaskTest.class);
    suite.addTestSuite(GraphMapParamsTest.class);
    suite.addTestSuite(GraphToPlotTest.class);
    suite.addTestSuite(GraphTraversionsTest.class);
    suite.addTestSuite(GraphSorterTest.class);
    suite.addTestSuite(HashMapDeBruijnGraphTest.class);
    suite.addTestSuite(KmerHashTest.class);
    suite.addTestSuite(KmerHashATest.class);
    suite.addTestSuite(AsyncKmerIterableTest.class);
    suite.addTestSuite(KmerIterableFactoryTest.class);
    suite.addTestSuite(KmerIteratorTest.class);
    suite.addTestSuite(LargeKDeBruijnGraphTest.class);
    suite.addTestSuite(LowKDeBruijnGraphTest.class);
    suite.addTestSuite(MergeNodesTest.class);
    suite.addTestSuite(NullTraversionsTest.class);
    suite.addTestSuite(PacBioCliTest.class);
    suite.addTestSuite(PacBioParamsTest.class);
    suite.addTestSuite(PacBioPathTest.class);
    suite.addTestSuite(PacBioStatisticsTest.class);
    suite.addTestSuite(PacBioTest.class);
    suite.addTestSuite(PairConstraintWriterTest.class);
    suite.addTestSuite(PairJoinerTest.class);
    suite.addTestSuite(PalindromeTrackerTest.class);
    suite.addTestSuite(PartialAlignmentTest.class);
    suite.addTestSuite(PathAlignerTest.class);
    suite.addTestSuite(PathTrackerTest.class);
    suite.addTestSuite(PreContigTest.class);
    suite.addTestSuite(ReadIteratorTest.class);
    suite.addTestSuite(ReadPairSourceTest.class);
    suite.addTestSuite(ReadPairSource454Test.class);
    suite.addTestSuite(ReadPairSourceMatePairTest.class);
    suite.addTestSuite(RecoverLinksTest.class);
    suite.addTestSuite(SequencesReaderIteratorTest.class);
    suite.addTestSuite(StringKmerTest.class);
    suite.addTestSuite(TraversionTest.class);

    suite.addTest(com.rtg.assembler.graph.AllTests.suite());

    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

