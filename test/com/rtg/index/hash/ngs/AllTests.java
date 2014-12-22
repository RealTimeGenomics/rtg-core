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
package com.rtg.index.hash.ngs;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory.
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.index.hash.ngs");

    suite.addTest(com.rtg.index.hash.ngs.instances.AllTests.suite());
    suite.addTest(com.rtg.index.hash.ngs.general.AllTests.suite());
    suite.addTest(com.rtg.index.hash.ngs.protein.AllTests.suite());

    suite.addTestSuite(ProteinTemplateHashLoopTest.class);
    suite.addTestSuite(ProteinIncrementalHashLoopTest.class);
    suite.addTestSuite(NgsHashLoopImplTest.class);
    suite.addTestSuite(ImplementHashFunctionTest.class);
    suite.addTestSuite(ImplementHashFunctionCloneTest.class);
    suite.addTestSuite(ReadCallImplementationTest.class);
    suite.addTestSuite(TemplateCallImplementationTest.class);
    suite.addTestSuite(PcrTemplateCallImplementationTest.class);
    suite.addTestSuite(FactoryUtilTest.class);
    suite.addTestSuite(MemScoreTest.class);
    suite.addTestSuite(ReadDecoderTest.class);
    suite.addTestSuite(ReadEncoderTest.class);
    suite.addTestSuite(ProteinNgsHashLoopTest.class);
    suite.addTestSuite(ProteinReadHashLoopTest.class);

    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}

