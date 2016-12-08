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
    //Diagnostic.developerLog("rlen: " + rlen + " start: " + zeroBasedStart);
//    System.err.println("read : " + Arrays.toString(read));
//    System.err.println("read len " + rlen);
//    System.err.println("template " + Arrays.toString(template));
//    System.err.println("start " + zeroBasedStart);
//    System.err.println("templateLength " + template.length);
//    System.err.println("templateEnd " + (zeroBasedStart + rlen));
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
