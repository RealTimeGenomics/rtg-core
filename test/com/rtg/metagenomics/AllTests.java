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
package com.rtg.metagenomics;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory. Run from the command
 * line with:<p>
 *
 * java -ea com.reeltwo.cartesian.AllTests
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.metagenomics");
    suite.addTest(com.rtg.metagenomics.matrix.AllTests.suite());
    suite.addTest(com.rtg.metagenomics.metasnp.AllTests.suite());
    suite.addTest(com.rtg.metagenomics.krona.AllTests.suite());
    suite.addTestSuite(BisectionSolverTest.class);
    suite.addTestSuite(BlockMappingTest.class);
    suite.addTestSuite(CompositionMetaPipelineCliTest.class);
    suite.addTestSuite(DummySolverTest.class);
    suite.addTestSuite(FragTest.class);
    suite.addTestSuite(FunctionalMetaPipelineCliTest.class);
    suite.addTestSuite(IdSetTest.class);
    suite.addTestSuite(LLineTest.class);
    suite.addTestSuite(LineSolverTest.class);
    suite.addTestSuite(LineTest.class);
    suite.addTestSuite(LinearInterpolationSolverTest.class);
    suite.addTestSuite(MetaPipelineParamsTest.class);
    suite.addTestSuite(MetagenomicsWrapperCliTest.class);
    suite.addTestSuite(MetagenomicsWrapperTaskTest.class);
    suite.addTestSuite(MinimizerTest.class);
    suite.addTestSuite(NewtonRaphsonSolverTest.class);
    suite.addTestSuite(SimpleTerminatorTest.class);
    suite.addTestSuite(SpeciesCliTest.class);
    suite.addTestSuite(SpeciesInfoTest.class);
    suite.addTestSuite(SpeciesLineDerivTest.class);
    suite.addTestSuite(SpeciesLineLinearDerivTest.class);
    suite.addTestSuite(SpeciesLineLinearTest.class);
    suite.addTestSuite(SpeciesLineTest.class);
    suite.addTestSuite(SpeciesMapTest.class);
    suite.addTestSuite(SpeciesParamsTest.class);
    suite.addTestSuite(SpeciesResultTest.class);
    suite.addTestSuite(SpeciesSdfSubsetTest.class);
    suite.addTestSuite(SpeciesStatisticsTest.class);
    suite.addTestSuite(SpeciesTargetCliTest.class);
    suite.addTestSuite(SpeciesTest.class);
    suite.addTestSuite(SubBlockResultTest.class);
    suite.addTestSuite(TargetTerminatorTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

