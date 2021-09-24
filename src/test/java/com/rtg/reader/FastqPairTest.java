/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.reader;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import com.rtg.mode.DnaUtils;

/**
 */
public class FastqPairTest {
  @Test
  public void leftRightAccessors() {
    final FastqPair fastqPair = new FastqPair(FastqSequenceTest.getFastq("name", "ACGT"), FastqSequenceTest.getFastq("name", "GGGGC"));
    assertEquals("ACGT", DnaUtils.bytesToSequenceIncCG(fastqPair.r1().getBases()));
    assertEquals("GGGGC", DnaUtils.bytesToSequenceIncCG(fastqPair.r2().getBases()));
  }
}