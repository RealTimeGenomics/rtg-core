/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    for (int i = 0; i < kmer.length(); ++i) {
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
    for (int i = 0; i < mKmerSize; ++i) {
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
