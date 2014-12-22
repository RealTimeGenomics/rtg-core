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

import com.rtg.util.LongUtils;


/**
 * Kmer specified by bits in a long.
 */
class KmerHash extends AbstractKmer {

  static long kmerToHashMin(Kmer kmer) {
    final Kmer rev = kmer.reverse();
    if (kmer.compareTo(rev) <= 0) {
      return kmerToHash(kmer);
    }
    return kmerToHash(rev);
  }

  static long kmerToHash(Kmer kmer) {
    assert kmer.length() <= 32;
    long hash = 0;
    for (int i = 0; i < kmer.length(); i++) {
      hash = hash << 2;
      hash |= kmer.nt(i) - 1;
    }
    return hash;
  }

  private final long mHash;
  private final int mKmerSize;
  private final long mMask;

  /**
   * @param hash the hash from the index.
   * @param size length of the kmer (in nt).
   */
  KmerHash(long hash, int size) {
    mHash = hash;
    mKmerSize = size;
    mMask = LongUtils.longMask(2 * size);
  }

  @Override
  public Kmer successor(byte nt) {
    final long sHash0 = (mHash << 2) | (nt - 1);
    return new KmerHash(sHash0 & mMask, mKmerSize);
  }

  @Override
  public Kmer predecessor(byte nt) {
    final int shift = 2 * mKmerSize - 2;
    final long nts = ((long) (nt - 1)) << shift;
    final long low = mHash >>> 2;
    final long pHash = low | nts;
    //System.err.println("mHash:" + Masks.bits(mHash, 64) + " nts:" + Masks.bits(nts, 64) + " low:" + Masks.bits(low, 64));
    return new KmerHash(pHash, mKmerSize);
  }

  @Override
  public Kmer minimalKmer() {
    final long reverse = reverseHash();
    //don't forget case when k == 32
    if (LongUtils.isLessThanUnsigned(reverse, mHash)) {
      return new KmerHash(reverse, mKmerSize);
    } else {
      return this;
    }
  }

  @Override
  public Kmer reverse() {
    final long reverse = reverseHash();
    return new KmerHash(reverse, mKmerSize);
  }

  private long reverseHash() {
    long reverse = 0;
    long forw = mHash;
    for (int i = 0; i < mKmerSize; i++) {
      final int nt = 3 - (int) (forw & 3);
      reverse = (reverse << 2) | nt;
      forw = forw >> 2;
    }
    return reverse;
  }

  @Override
  public int length() {
    return mKmerSize;
  }

  @Override
  public byte nt(int pos) {
    final int ix = 2 * (mKmerSize - pos - 1);
    final long b = (mHash >> ix) & 3;
    return (byte) (b + 1);
  }
}
