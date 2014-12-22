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
package com.rtg.util.array.intindex;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this directory. Run from the command
 * line with:<p>
 *
 * java -ea com.reeltwo.util.array.big.AllTests
 *
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.util.array.intindex");
    suite.addTestSuite(IntArrayTest.class);
    suite.addTestSuite(IntChunksTest.class);
    suite.addTestSuite(IntCreateTest.class);
    suite.addTestSuite(IntIndexTest.class);
    suite.addTestSuite(SmallChunksTest.class);
    return suite;
  }
}
