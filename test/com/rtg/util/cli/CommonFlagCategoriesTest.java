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
package com.rtg.util.cli;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;


/**
 */
public class CommonFlagCategoriesTest extends TestCase {

  public void testSetCategories() {
    CFlags flags = new CFlags();
    flags.registerOptional("a", "a").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerOptional("b", "b").setCategory(CommonFlagCategories.REPORTING);
    flags.registerOptional("c", "c").setCategory(CommonFlagCategories.SENSITIVITY_TUNING);
    flags.registerOptional("d", "d").setCategory(CommonFlagCategories.UTILITY);
    flags.registerOptional("e", "e").setCategory(CommonFlagCategories.FILTERING);
    CommonFlagCategories.setCategories(flags);
    assertEquals("Utility", flags.getFlag("help").getCategory());
    TestUtils.containsAll(flags.getUsageString(), "File Input/Output", "Reporting", "Sensitivity Tuning", "Utility");
  }
}
