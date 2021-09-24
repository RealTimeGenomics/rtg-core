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
import com.rtg.util.integrity.Exam;

/**
 * Kmer specified by bits in multiple longs, handles k &gt; 32.
 * <br>
 * The nucleotides are laid out as follows in the array of longs.
 * This example shows the positions of the nucleotides which are counted from 0.
 * The kmer length is 41 which means 82 bits are used, 64 in long[0] and 18 in long[1]
 * <pre>
 * long[0]                                          long[1]
 * +-----------------------------------------+      +-----------------------------------------+
 *  18 19 20 ........................39 40 41        -----------------------------0 1 2 ... 17
 * </pre>
 * The high order bits in long[1] are zero.<br>
 * This can handle higher values of k which may require more longs.
 */
class KmerHashA extends AbstractKmer {

  private static final int NT_IN_LONG = Long.SIZE / 2;

  static long[] kmerToHashMin(Kmer kmer) {
    final Kmer rev = kmer.reverse();
    if (kmer.compareTo(rev) <= 0) {
      return kmerToHash(kmer);
    }
    return kmerToHash(rev);
  }

  static long[] kmerToHash(Kmer kmer) {
    final int nts = kmer.length();
    final int length = (nts  + NT_IN_LONG - 1) / NT_IN_LONG;
    final int excess = nts % NT_IN_LONG;
    final long[] hash = new long[length];
    for (int i = 0, j = length - 1; i < kmer.length();) {
      hash[j] = (hash[j] << 2) | (kmer.nt(i) - 1);
      ++i;
      if (i % NT_IN_LONG == excess) {
        --j;
      }
    }
    return hash;
  }

  //  static long[] min(final long[] a, final long[] b) {
  //    assert a.length == b.length;
  //    for (int i = 0; i < a.length; ++i) {
  //      final long la = a[i];
  //      final long lb = b[i];
  //      if (la < lb) {
  //        return a;
  //      }
  //      if (la > lb) {
  //        return b;
  //      }
  //    }
  //    return a; //a equals b
  //  }

  static int lastBits(int bits) {
    return (bits - 1) % Long.SIZE + 1;
  }

  private final long[] mHash;
  private final int mKmerSize;
  private final long mLastMask;
  private final int mLastBits;

  /**
   * @param hash the hash from the index.
   * @param size length of the kmer (in nt).
   */
  KmerHashA(long[] hash, int size) {
    mHash = hash;
    mKmerSize = size;
    mLastBits = lastBits(2 * size);
    mLastMask = LongUtils.longMask(mLastBits);
  }

  @Override
  public Kmer successor(final byte nt) {
    final int length = mHash.length;
    final long[] succ = new long[length];
    int carry = nt - 1;
    for (int i = 0; i < length; ++i) {
      final long l = mHash[i];
      final int next = (int) (l >>> (Long.SIZE - 2));
      succ[i] = (l << 2) | carry;
      carry = next;
    }
    succ[length - 1] &= mLastMask;
    return new KmerHashA(succ, mKmerSize);
  }

  @Override
  public Kmer predecessor(byte nt) {
    final int length = mHash.length;
    //System.err.println(length);
    final long[] pred = new long[length];
    final long l0 = mHash[length - 1];
    long carry = l0;
    final long hi = ((long) (nt - 1)) << (mLastBits - 2);
    final long lo = l0 >>> 2;
    pred[length - 1] = hi | lo;
    for (int i = length - 2; i >= 0; --i) {
      final long l = mHash[i];
      final long ca = carry << (Long.SIZE - 2); //because of this shift we dont need to mask carry
      pred[i] = ca | l >>> 2;
      carry = l;
    }
    return new KmerHashA(pred, mKmerSize);
  }

  static boolean isLessThanUnsigned(long[] n1, long[] n2) {
    assert n1.length == n2.length;
    for (int k = 0; k < n1.length; ++k) {
      if (n1[k] != n2[k]) {
        return LongUtils.isLessThanUnsigned(n1[k], n2[k]);
      }
    }
    return false;
  }

  @Override
  public Kmer minimalKmer() {
    final long[] reverse = reverseHash();
    //don't forget case when k == 32
    if (isLessThanUnsigned(reverse, mHash)) {
      return new KmerHashA(reverse, mKmerSize);
    } else {
      return this;
    }
  }

  @Override
  public Kmer reverse() {
    final long[] reverse = reverseHash();
    return new KmerHashA(reverse, mKmerSize);
  }

  private long[] reverseHash() {
    if (mHash.length == 0) {
      return mHash;
    }
    final long[] reverse = new long[mHash.length];
    final int excess = mKmerSize % NT_IN_LONG;
    int i = 0;
    int j = 0;
    int k = reverse.length - 1;
    long f = mHash[j];
    while (true) {
      final int nt = 3 - (int) (f & 3);
      reverse[k] = (reverse[k] << 2) | nt;
      f = f >> 2;
      if (++i >= mKmerSize) {
        break;
      }
      if (i % NT_IN_LONG == 0) {
        f = mHash[++j];
      }
      if (i % NT_IN_LONG == excess) {
        --k;
      }
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
    final long b = (mHash[ix >> 6] >> ix) & 3;
    return (byte) (b + 1);
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertTrue(mKmerSize > 32);
    Exam.assertTrue(mHash.length * 32 >= mKmerSize);
    return true;
  }
}
