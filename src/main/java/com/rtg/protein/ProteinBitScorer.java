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
package com.rtg.protein;


/**
 * Takes read and template as protein and works out a rough lower bound
 * of the percent identity.
 *
 */
final class ProteinBitScorer {

  /** Maximum length insert/delete that will still pass the <code>%identity.</code> */
  private final int mMaxIndel;

  private static final int NUM_BITS = 5;
  private final long[] mTemplate = new long[NUM_BITS];
  private final long[] mRead = new long[NUM_BITS];

  private byte[] mPrevTemplate = null;
  private int mPrevZeroBasedStart = Integer.MAX_VALUE;

  /**
   * Contains number of equal bits for shifts of <code>-mMaxIndel..+mMaxIndel</code>.
   * The last entry is the number of equal bits in the bitwise OR of all shifts.
   */
  private final int[] mResults;

  /**
   * Create a protein percent identity estimator.
   * @param maxIndel maximum length insert/delete that will have scores calculated
   */
  ProteinBitScorer(final int maxIndel) {
    mMaxIndel = maxIndel;
    mResults = new int[2 * mMaxIndel + 2];
  }

  public int getMaxIndelLength() {
    return mMaxIndel;
  }

  /**
   * Estimate the number of identical amino acids.
   *
   * @param read amino acids
   * @param zeroBasedReadStart start position in read to process
   * @param rlen amount of read to process
   * @param template template amino acids
   * @param zeroBasedStart start position within template
   * @return an array of bit-scores, for shifts of <code>-n .. +n</code>, followed
   *     by the bit-score of the OR of all the shifts.
   */
  public int[] calculateFastScore(final byte[] read, int zeroBasedReadStart, final int rlen, final byte[] template, final int zeroBasedStart) {
//    System.err.println("read: " + Protein.bytesToProtein(read, 0, read.length));
//    System.err.println("read [" + zeroBasedReadStart + "," + (zeroBasedReadStart + rlen) + "):" + Protein.bytesToProtein(read, zeroBasedReadStart, rlen));
//    System.err.println("template " + Protein.bytesToProtein(template, 0, template.length));
//    System.err.println("template [" + zeroBasedStart + "," + (zeroBasedStart + rlen) + "): " + Protein.bytesToProtein(template, zeroBasedStart, rlen));
    fillup(read, zeroBasedReadStart, rlen, 0, mRead);
    // we include an extra mMaxIndel amino acids on each side of the template.
    if (zeroBasedStart != mPrevZeroBasedStart || template != mPrevTemplate || rlen != 0) {
      fillup(template, zeroBasedStart - mMaxIndel, rlen + 2 * mMaxIndel, 127, mTemplate);
      mPrevZeroBasedStart = zeroBasedStart;
      mPrevTemplate = template;
    }

    long finalBits = 0;
    for (int shift = 0; shift <= 2 * mMaxIndel; ++shift) {
      long eqBits = (1L << rlen) - 1;
      for (int bit = 0; bit < NUM_BITS; ++bit) {
        final long thisEq = ~(mTemplate[bit] ^ (mRead[bit] << shift));
        eqBits &= thisEq >> shift;
      }
      mResults[shift] = Long.bitCount(eqBits);
      finalBits |= eqBits;
    }
    mResults[mResults.length - 1] = Long.bitCount(finalBits);
    return mResults;
  }

  /**
   * Extracts bit <code>bitNum</code> from each byte in
   * <code>data[start..start+len-1]</code> into <code>store[bitNum]</code>.
   *
   * @param data an array of amino acids.
   * @param start start position (can be negative).
   * @param len number of bytes to process (can go past end of data).
   * @param unknown an impossible amino acid.  This is used outside <code>data</code>.
   * @param store where the result bit vectors are returned.
   */
  void fillup(final byte[] data, final int start, final int len, final int unknown, final long[] store) {
    long a = 0, b = 0, c = 0, d = 0, e = 0;
    final int end = start + len;
    for (int i = start; i < end; ++i) {
      final byte aa = 0 <= i && i < data.length ? data[i] : (byte) unknown;
      a = a << 1 | (aa & 1);
      b = b << 1 | ((aa >> 1) & 1);
      c = c << 1 | ((aa >> 2) & 1);
      d = d << 1 | ((aa >> 3) & 1);
      e = e << 1 | ((aa >> 4) & 1);
    }
    store[0] = a;
    store[1] = b;
    store[2] = c;
    store[3] = d;
    store[4] = e;
  }
}
