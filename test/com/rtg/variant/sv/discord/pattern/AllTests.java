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

package com.rtg.variant.sv.discord.pattern;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 *         Date: 16/03/12
 *         Time: 11:08 AM
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant.sv.discord.pattern");

    suite.addTestSuite(BreakpointAltTest.class);
    suite.addTestSuite(BreakpointPatternParamsTest.class);
    suite.addTestSuite(BreakpointStoreTest.class);
    suite.addTestSuite(DeletionOverlapFilterTest.class);
    suite.addTestSuite(SvPatternsCliTest.class);
    suite.addTestSuite(SvPatternsTaskTest.class);
    suite.addTestSuite(VcfBreakpointTest.class);

    return suite;
  }
}
