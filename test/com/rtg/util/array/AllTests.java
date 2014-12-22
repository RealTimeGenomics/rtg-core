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
package com.rtg.util.array;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory. Run from the command
 * line with:<p>
 *
 * java -ea com.reeltwo.util.array.AllTests
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.util.array");
    suite.addTestSuite(ArrayHandleTest.class);
    suite.addTestSuite(ArrayTypeTest.class);
    suite.addTestSuite(DummyIndexTest.class);
    suite.addTestSuite(IndexTypeTest.class);
    suite.addTestSuite(WrappedIntArrayTest.class);
    suite.addTestSuite(ArrayUtilsTest.class);
    suite.addTestSuite(SingleValueIntArrayTest.class);

    suite.addTest(com.rtg.util.array.bitindex.AllTests.suite());
    suite.addTest(com.rtg.util.array.byteindex.AllTests.suite());
    suite.addTest(com.rtg.util.array.intindex.AllTests.suite());
    suite.addTest(com.rtg.util.array.longindex.AllTests.suite());
    suite.addTest(com.rtg.util.array.objectindex.AllTests.suite());
    suite.addTest(com.rtg.util.array.shortindex.AllTests.suite());
    suite.addTest(com.rtg.util.array.packedindex.AllTests.suite());
    suite.addTest(com.rtg.util.array.zeroindex.AllTests.suite());
    return suite;
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }
}
