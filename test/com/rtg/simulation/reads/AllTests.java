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
package com.rtg.simulation.reads;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.simulation.reads");
    suite.addTestSuite(CgSimCliTest.class);
    suite.addTestSuite(CompleteGenomicsMachineTest.class);
    suite.addTestSuite(DummyIlluminaMachineTest.class);
    suite.addTestSuite(DummyMachineTest.class);
    suite.addTestSuite(ErrorMachineTest.class);
    suite.addTestSuite(FastaReadWriterTest.class);
    suite.addTestSuite(FastqReadWriterTest.class);
    suite.addTestSuite(FilteringFragmenterTest.class);
    suite.addTestSuite(FourFiveFourPairedEndMachineTest.class);
    suite.addTestSuite(FourFiveFourSingleEndMachineTest.class);
    suite.addTestSuite(GenomeFragmenterTest.class);
    suite.addTestSuite(IlluminaPairedEndMachineTest.class);
    suite.addTestSuite(IlluminaSingleEndMachineTest.class);
    suite.addTestSuite(IonTorrentSingleEndMachineTest.class);
    suite.addTestSuite(ReadSimCliTest.class);
    suite.addTestSuite(ReadSimValidatorTest.class);
    suite.addTestSuite(SdfReadWriterTest.class);
    suite.addTestSuite(TaxonomyDistributionTest.class);
    suite.addTestSuite(UnknownBaseReadWriterTest.class);

    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
