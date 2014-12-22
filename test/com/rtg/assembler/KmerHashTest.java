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
public class KmerHashTest extends AbstractKmerTest {

  @Override
  protected Kmer kmer(String str) {
    final int size = str.length();
    final Kmer k0 = new StringKmer(str);
    final long hash = KmerHash.kmerToHash(k0);
    return new KmerHash(hash, size);
  }

  public void testKmerToHashMin() {
    assertEquals(26, KmerHash.kmerToHash(kmer("ACGG")));
    assertEquals(26, KmerHash.kmerToHashMin(kmer("ACGG")));
    assertEquals(91, KmerHash.kmerToHash(kmer("CCGT")));
    assertEquals(26, KmerHash.kmerToHashMin(kmer("CCGT")));
  }
}
