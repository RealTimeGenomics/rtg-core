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
package com.rtg.simulation.snpsim;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class GenDiffModeTest extends TestCase {

  public void test() {
    TestUtils.testEnum(GenDiffMode.class, "[BOTH_SAME, DIFFERENT, FIRST_ONLY, TWIN_ONLY]");
  }

  public void testDiffModeArray() {
    final String[] gen = {"ac", "ct", "gt"};
    TestUtils.assertEquals(new String[] {"ac", "ct"}, GenDiffMode.diffModeArray(gen, GenDiffMode.BOTH_SAME));
    TestUtils.assertEquals(new String[] {"ac", "ct", "gt"}, GenDiffMode.diffModeArray(gen, GenDiffMode.DIFFERENT));
    TestUtils.assertEquals(new String[] {"ac", "ct", "ac"}, GenDiffMode.diffModeArray(gen, GenDiffMode.FIRST_ONLY));
    TestUtils.assertEquals(new String[] {"ac", "ac", "gt"}, GenDiffMode.diffModeArray(gen, GenDiffMode.TWIN_ONLY));
  }
}
