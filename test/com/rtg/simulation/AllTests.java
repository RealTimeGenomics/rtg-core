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
package com.rtg.simulation;


import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * AllTests for reader package
 *
 */
public class AllTests extends TestCase {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.simulation");
    suite.addTest(com.rtg.simulation.genome.AllTests.suite());
    suite.addTest(com.rtg.simulation.reads.AllTests.suite());
    suite.addTest(com.rtg.simulation.variants.AllTests.suite());
    suite.addTest(com.rtg.simulation.snpsim.AllTests.suite());
    suite.addTest(com.rtg.simulation.sv.AllTests.suite());
    suite.addTest(com.rtg.simulation.cnv.AllTests.suite());
    suite.addTestSuite(DwgsimReadNameParserTest.class);
    suite.addTestSuite(MutatedOffsetsTest.class);
    suite.addTestSuite(MutatedReferenceReadNameParserTest.class);
    suite.addTestSuite(NewReadNameParserTest.class);
    suite.addTestSuite(NewestReadNameParserTest.class);
    suite.addTestSuite(OldReadNameParserTest.class);
    suite.addTestSuite(ReadMappingAccuracyParamsTest.class);
    suite.addTestSuite(ReadMappingAccuracyReadStatsTest.class);
    suite.addTestSuite(ReadMappingAccuracyTest.class);
    suite.addTestSuite(ReadMappingRocTest.class);
    suite.addTestSuite(SimulationUtilsTest.class);
    suite.addTestSuite(SimulatedReadNameParserFactoryTest.class);
    suite.addTestSuite(SoftClipCigarParserTest.class);
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

