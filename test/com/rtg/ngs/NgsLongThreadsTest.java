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
package com.rtg.ngs;

import com.rtg.ngs.NgsTestUtils.NgsFilterPartlyParams;
import com.rtg.ngs.NgsTestUtils.TestParams;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test threading for long reads
 */
public class NgsLongThreadsTest extends NgsLongTest {

  @Override
  void checkLong(NgsMaskParams mask, TestParams test, NgsFilterPartlyParams filterParams, boolean containsOnly) throws Exception {
    NgsTestUtils.checkLongThreads(mask, test, filterParams, containsOnly);
  }

  public static Test suite() {
    TestSuite suite = new TestSuite();

    suite.addTestSuite(NgsLongThreadsTest.class);
    return suite;
  }
}
