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

package com.rtg.assembler;

/**
 */
public class HashMapDeBruijnGraphTest extends AbstractDeBruijnGraphTest {

  @Override
  HashMapDeBruijnGraph getDeBruijnGraph(final KmerMockFactory kit, long size, int kMerSize) {
    return new HashMapDeBruijnGraph(kit);
  }

  public void testAddContains() {
    aTest(4, 5, "AAAAA", "AAACC", "ATACC");
  }

  public void testAddContainsBig() {
    final int kmerSize = 35;
    final String sNot = big("AAAAA", kmerSize);
    final String s0 = big("AAACC", kmerSize);
    final String s1 = big("ATACC", kmerSize);
    aTest(64, kmerSize, sNot, s0, s1);
  }
}
