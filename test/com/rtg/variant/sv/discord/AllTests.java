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
package com.rtg.variant.sv.discord;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * Test class for all tests in this package.
 */
public class AllTests extends TestSuite {

  public static Test suite() {
    final TestSuite suite = new TestSuite("com.rtg.variant.sv");
    suite.addTest(BreakpointConstraintTest.suite());
    suite.addTestSuite(BreakpointGeometryTest.class);
    suite.addTestSuite(BreakpointPositionTest.class);
    suite.addTestSuite(BedDiscordantOutputFormatterTest.class);
    suite.addTestSuite(DebugComparatorTest.class);
    suite.addTestSuite(DebugDiscordantOutputFormatterTest.class);
    suite.addTestSuite(DiscordantReadSetTest.class);
    suite.addTestSuite(DiscordantToolTest.class);
    suite.addTestSuite(DiscordantToolCliTest.class);
    suite.addTestSuite(DiscordantToolParamsTest.class);
    suite.addTestSuite(DiscordantToolStatisticsTest.class);
    suite.addTestSuite(DummyBreakpointGeometryTest.class);
    suite.addTestSuite(FlipTest.class);
    suite.addTestSuite(IntervalTest.class);
    suite.addTestSuite(OrientationTest.class);
    suite.addTestSuite(ReorderingDebugOutputTest.class);
    suite.addTestSuite(VcfDiscordantOutputFormatterTest.class);
    suite.addTestSuite(SmartBedWriterTest.class);
    suite.addTestSuite(SmartVcfWriterTest.class);
    suite.addTest(com.rtg.variant.sv.discord.pattern.AllTests.suite());
    return suite;
  }
}
