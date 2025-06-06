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

import java.util.Arrays;

import com.rtg.mode.DnaUtils;

/**
 */
public class KmerHashATest extends AbstractKmerTest {

  @Override
  protected Kmer kmer(String str) {
    final int size = str.length();
    final Kmer k0 = new StringKmer(str);
    final long[] hash = KmerHashA.kmerToHash(k0);
    return new KmerHashA(hash, size);
  }

  //                                     0         1         2         3         4         5         6         7
  //                                     01234567890123456789012345678901234567890123456789012345678901234567890123456789
  private static final String TEST_NT = "TGCACGTGTACGTCACCGTGGTTTGGCCCGCCACCACTGGCCCCCGCCGGACCAGTGTCGTGCACGTGCAAACGGTGTCAA";

  public void testKmerToHash() {
    final String nt0 = TEST_NT;
    for (int j = 0; j < nt0.length(); ++j) {
      final String nt = nt0.substring(0, j);
      final Kmer k0 = new StringKmer(nt);
      final long[] hash = KmerHashA.kmerToHash(k0);
      final Kmer k = new KmerHashA(hash, nt.length());
      for (int i = 0; i < nt.length(); ++i) {
        assertEquals("j=" + j + " i=" + i + " : " + Arrays.toString(hash), nt.charAt(i), DnaUtils.getBase(k.nt(i)));
      }
    }
  }

  public void testReverseHash() {
    final String nt0 = TEST_NT;
    for (int j = 0; j < nt0.length(); ++j) {
      final String nt = nt0.substring(0, j);
      final Kmer k0 = new StringKmer(nt);
      final long[] hash = KmerHashA.kmerToHash(k0);
      final Kmer k = new KmerHashA(hash, nt.length());
      assertEquals("j=" + j, k, k.reverse().reverse());
      final String rc = DnaUtils.reverseComplement(nt);
      final Kmer rck = new KmerHashA(KmerHashA.kmerToHash(new StringKmer(rc)), rc.length());
      assertEquals(rck, k.reverse());
    }
  }

  public void testLastBits() {
    assertEquals(0, KmerHashA.lastBits(0));
    assertEquals(1, KmerHashA.lastBits(1));
    assertEquals(64, KmerHashA.lastBits(64));
    assertEquals(1, KmerHashA.lastBits(65));
    assertEquals(64, KmerHashA.lastBits(128));
    assertEquals(1, KmerHashA.lastBits(129));
  }

  //  public void testMin() {
  //    checkMin(new long[] {}, new long[] {});
  //    checkMin(new long[] {0}, new long[] {0});
  //    checkMin(new long[] {0}, new long[] {1});
  //    checkMin(new long[] {0, 1}, new long[] {0, 1});
  //    checkMin(new long[] {0, 1}, new long[] {0, 2});
  //    checkMin(new long[] {0, 1}, new long[] {1, 1});
  //  }
  //
  //  private void checkMin(long[] la, long[] lb) {
  //    assertTrue(la == KmerHashA.min(la, lb));
  //    assertTrue(la == KmerHashA.min(la, la));
  //    assertTrue(lb == KmerHashA.min(lb, lb));
  //    if (!Arrays.equals(la, lb)) {
  //      assertTrue(la == KmerHashA.min(lb, la));
  //    }
  //  }

  private static final String LONG_KMER;
  static {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < 18; ++i) {
      sb.append("GACT");
    }
    LONG_KMER = sb.toString();
  }
  private static final String VERY_LONG_KMER = LONG_KMER + LONG_KMER;

  public void testKmerToHashMin() {
    checkKmerToHashMin("ACTG");
    checkKmerToHashMin(LONG_KMER);
  }

  private void checkKmerToHashMin(final String str) {
    final Kmer kmer = kmer(str);
    final long[] l0 = KmerHashA.kmerToHashMin(kmer);
    final long[] lf = KmerHashA.kmerToHash(kmer);
    final long[] lr = KmerHashA.kmerToHash(kmer.reverse());
    assertTrue(Arrays.equals(l0, lf) || Arrays.equals(l0, lr));
  }

  @Override
  public void testSuccessor() {
    super.testSuccessor();
    checkSuccessor(LONG_KMER);
    checkSuccessor(VERY_LONG_KMER);
  }

  @Override
  public void testPredecessor() {
    super.testPredecessor();
    checkPredecessor(LONG_KMER);
    checkPredecessor(VERY_LONG_KMER);
  }

  @Override
  public void testMinimal() {
    super.testMinimal();
    checkMinimal(LONG_KMER);
    checkMinimal(VERY_LONG_KMER);
  }

  public void testIsLessThanUnsigned() {
    checkIsLessThanUnsigned(new long[] {0}, new long[] {1});
    checkIsLessThanUnsigned(new long[] {0, 1}, new long[] {0, 2});
    checkIsLessThanUnsigned(new long[] {0, 1}, new long[] {1, 1});
  }

  private void checkIsLessThanUnsigned(long[] la, long[] lb) {
    assertTrue(KmerHashA.isLessThanUnsigned(la, lb));
    assertFalse(KmerHashA.isLessThanUnsigned(lb, la));
    assertFalse(KmerHashA.isLessThanUnsigned(la, la));
    assertFalse(KmerHashA.isLessThanUnsigned(lb, lb));
  }

  public void testReverseBug() {
    assertEquals("GGTTT", kmer("AAACC").reverse().toString());
  }
}
