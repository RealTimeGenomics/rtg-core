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

package com.rtg.metagenomics;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class SpeciesStatisticsTest extends TestCase {

  public void testStats() {
    SpeciesStatistics s = new SpeciesStatistics(null);
    s.mRichness = 5;
    s.mShannon = 2.5;
    s.mPielou = 1;
    s.mInvSimpson = 0.3;
    String res = s.getStatistics();
    TestUtils.containsAll(res, "Richness", "5", "Shannon", "2.5", "Pielou", "1", "Inverse Simpson", "0.3");
  }
}
